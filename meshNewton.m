function [z,triEn, allStats] = meshNewton(xref, t, P2PVtxIds, P2PDst, y, nIter, p2p_weight, energy_type, energy_param, hession_proj)

fR2C = @(x) complex(x(:,1), x(:,2));
fC2R = @(x) [real(x) imag(x)];
if ~exist('splsolver_imp','file')
    addpath( genpath([pwd '\splsolver']) );
end

r=0.5;
barrier_param=0.00001;
ub=0.45; db=-0.45;
lb=-0.45; rb=0.45;

scal=1;tmp=0.00001;
fB=@(x,param_e) 1/(scal*x+sqrt((scal*x).^2+param_e))*tmp;
dfB=@(x,param_e) -scal./((scal*x).*(sqrt((scal*x).^2+param_e)+(scal*x))+param_e)*tmp;
ddfB=@(x,param_e) scal^2./sqrt(((scal*x).^2+param_e).^3)*tmp;
eB=0.00001;

%barrier_param=0;

if ~exist('energy_type', 'var'), hession_proj = 'SymmDirichlet'; end
if ~exist('hession_proj', 'var'), hession_proj = 'KP'; end
if ~exist('energy_param', 'var'), energy_param = 1; end


if isreal(xref), xref = fR2C(xref); end
if isreal(y), y = fR2C(y); end
if ~isempty(P2PDst)
    if isreal(P2PDst), P2PDst = fR2C(P2PDst); end
end
nf = size(t, 1);
nv = size(y, 1);
% 
% if isempty(P2PVtxIds)
%     [~,P2PVtxIds] = min(abs(y-mean(y)));
%     P2PDst = y(P2PVtxIds);
% end

%%
if size(xref,1)==nv
    xref = xref(t(:,2:3)) - xref(t(:,1)); % implicit the x1 = 0, for each triangle
end


%% initilization
isometric_energyies = [ "SymmDirichlet", "ExpSD", "AMIPS", "SARAP", "HOOK", "ARAP", "BARAP", "BCONF"];
hessian_projections = [ "NP", "KP", "FP4", "FP6", "CM" ];



findStringC = @(s, names) find(strcmpi(s, names), 1) - 1;
mexoption = struct('energy_type', findStringC(energy_type, isometric_energyies), ...
                   'hessian_projection', findStringC(hession_proj, hessian_projections), ...
                   'energy_param', energy_param, 'verbose', 0);

z = y;

Areas = (real(xref(:,1)).*imag(xref(:,2)) - imag(xref(:,1)).*real(xref(:,2)))/2;
D2 = -1i/4*(xref*[1 0 -1; -1 1 0])./Areas;
D2t = D2.';
D = sparse(repmat(1:nf,3,1)', t, D2);


% fIsoEnergyFromFzGz = @(fz, gz) meshIsometricEnergy(D2, fz, gz, Areas, energy_type, energy_param);
fIsoEnergyFromFzGz = @(fz, gz) meshIsometricEnergyC(fz, gz, D2t, Areas, mexoption);

%fDeformEnergy = @(z) fIsoEnergyFromFzGz(conj(D*conj(z)), D*z) + p2p_weight*norm(z(P2PVtxIds)-P2PDst)^2;
%fDeformEnergy = @(z) fIsoEnergyFromFzGz(conj(D*conj(z)), D*z) + p2p_weight*norm(z(P2PVtxIds)-P2PDst)^2+barrier_param*sum(1./((r^2-real(z).^2))+1./((r^2-imag(z).^2)));
fDeformEnergy = @(z) fIsoEnergyFromFzGz(conj(D*conj(z)), D*z) + p2p_weight*norm(z(P2PVtxIds)-P2PDst)^2+2*sum(fB(real(z)-lb,eB)+fB(rb-real(z),eB)+fB(imag(z)-db,eB)+fB(ub-imag(z),eB));
%% initialization, get sparse matrix pattern
[xmesh, ymesh] = meshgrid(1:6, 1:6);
t2 = [t t+nv]';
Mi = t2(xmesh(:), :);
Mj = t2(ymesh(:), :);

 H = sparse(Mi, Mj, 1, nv*2, nv*2);  % only pattern is needed
%L = laplacian(y, t, 'uniform');
%L=cotLaplace(fC2R(y),t);
%H = [L L; L L];

% nonzero indices of the matrix
Hnonzeros0 = zeros(nnz(H),1);
Hnonzeros2=zeros(nnz(H),1);
P2PxId = uint64( [P2PVtxIds P2PVtxIds+nv] );
idxDiagH = ij2nzIdxs(H, P2PxId, P2PxId);
idxDiagH2 = ij2nzIdxs(H,uint64(1:nv*2),uint64(1:nv*2));
Hnonzeros0(idxDiagH) = p2p_weight*2;
nzidx = ij2nzIdxs(H, uint64(Mi), uint64(Mj));

G0 = zeros(nv*2,1);

solver = splsolver(H, 'ldlt');

allStats = zeros(nIter+1, 8); % statistics
allStats(1, [5 7 8]) = [0 norm(z(P2PVtxIds)-P2PDst)^2 fDeformEnergy(z)];

%% main loop
g2GIdx = uint64(t2);
for it=1:nIter
    tt = tic;

    fz = conj(D*conj(z)); % equivalent but faster than conj(D)*z;
    gz = D*z;
    
    assert( all( abs(fz)>abs(gz) ) );

    [en, g, hs] = fIsoEnergyFromFzGz(fz, gz);
    enIso=en;
    en = en + p2p_weight*norm(z(P2PVtxIds)-P2PDst)^2;
    dp2p = z(P2PVtxIds) - P2PDst;
    G0(P2PxId) = [real(dp2p); imag(dp2p)]*p2p_weight*2;
    G = myaccumarray(g2GIdx, g, G0);
    Hnonzeros = myaccumarray( nzidx, hs, Hnonzeros0 );
    
    Rz=real(z); Iz=imag(z);
%    G1=barrier_param*[2*Rz./(r^2-Rz.*Rz).^2;2*Iz./(r^2-Iz.*Iz).^2];
    G1=2*[-dfB(r-Rz,eB)+dfB(Rz+r,eB);-dfB(r-Iz,eB)+dfB(Iz+r,eB)];
    G=G+G1;
%    Hnonzeros2(idxDiagH2)=barrier_param*[2*(r^2+3*Rz.*Rz)./(r^2-Rz.*Rz).^3;2*(r^2+3*Iz.*Iz)./(r^2-Iz.*Iz).^3];
    Hnonzeros2(idxDiagH2)=2*[ddfB(r-Rz,eB)+ddfB(Rz+r,eB);ddfB(r-Rz,eB)+ddfB(Rz+r,eB)];
    Hnonzeros=Hnonzeros+Hnonzeros2;
%    en=en+barrier_param*sum(1./((r^2-real(z).^2))+1./((r^2-imag(z).^2)));
    en=en+2*sum(fB(real(z)-lb,eB)+fB(rb-real(z),eB)+fB(imag(z)-db,eB)+fB(ub-imag(z),eB));
    %% Newton
    dz = solver.refactor_solve(Hnonzeros, -G);
    dz = fR2C( reshape(dz, [], 2) );
    
    %% orientation preservation
    ls_t = min( maxtForPositiveArea2( fz, gz, conj(D*conj(dz)), D*dz )*0.9, 1 );
    
    %% line search energy decreasing
    fMyFun = @(t) fDeformEnergy( dz*t + z );
    normdz = norm(dz);

    dgdotfz = dot( G, [real(dz); imag(dz)] );
    
    ls_alpha = 0.2; ls_beta = 0.5;
    fQPEstim = @(t) en+ls_alpha*t*dgdotfz;

    e_new = fMyFun(ls_t);
    while ls_t*normdz>1e-12 && e_new > fQPEstim(ls_t)
        ls_t = ls_t*ls_beta;
        e_new = fMyFun(ls_t);
    end
    en = e_new;
    
    fprintf('\nit: %3d, en: %.3e, runtime: %.3fs, ls: %.2e', it, en, toc(tt), ls_t*normdz);
    assert( all( signedAreas(dz*ls_t + z, t)>0 ) );
    
    %% update
    z = dz*ls_t + z;
    
    %% stats
    allStats(it+1, [5 7 8]) = [toc(tt)*1000 norm(z(P2PVtxIds)-P2PDst)^2 en];
    
    
    if norm(G)<1e-4, break; end
    if norm(dz)*ls_t<1e-10, break; end
end

fz = conj(D*conj(z)); % equivalent but faster than conj(D)*z;
gz = D*z;
triEn = energyForEverySingleTriangle(fz,gz,Areas);

allStats = allStats(1:it+1, :);
allStats(:,7:8) = allStats(:,7:8)/sum(Areas);

