%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef OptimProblemIsoDist < OptimProblem
    properties
        % mesh (keep here?)
        V;
        F;
        dim;
        n_vert;
        n_tri;
        areas;
        
        % internal variables
        Tx;
        R;
        f_val;
        Tx_grad;
        flips;
        localHess;
        n_flip;
        
        % parameters
        maxStepTol = 1e-8;
        initArapIter;
        
        % bd projection parameters (for initialization)
        proj_bd_K = 10; % bounded distortion constant
        proj_bd_lb = -1; % lower bound on SVs (-1 = disabled)
        proj_bd_ub = -1; % upper bound on SVs (-1 = disabled)
        proj_bd_iter_max = 1000; % maximal number of BD projection iterations
        proj_bd_tol_err = 1e-10; % tolerance for stopping BD projection iterations
        proj_bd_verbose = true; % flag to showing iteration statistics
    end
    
    methods
        function obj = OptimProblemIsoDist(V, F, eq_lhs, eq_rhs, V0, initArapIter)
            % copy
            obj.V = V;
            obj.F = F;
            obj.eq_lhs = eq_lhs;
            obj.eq_rhs = eq_rhs;
            obj.dim = size(F,2)-1;
            obj.n_vert = size(V,1);
            obj.n_tri = size(F,1);
            if exist('initArapIter','var')
                obj.initArapIter = initArapIter;
            else
                obj.initArapIter = 1;
            end
            % report
            obj.report(1,'Constructing %s (dim: %d   #vert: %d   #elem: %d)\n', class(obj), obj.dim, obj.n_vert, obj.n_tri);
            % compute transformations
            [obj.T,obj.areas] = computeMeshTranformationCoeffsMex(F, V);
            % set initial configuration
            obj.x0 = obj.initVertices(V0);
            % set quadratic proxy
            obj.H = obj.setQuadraticProxy();
            % finish construction
            obj.initProblem();
        end
        
        function x0 = initVertices(obj, V0)
            % feasible initialization -- provide input or a single local-global step
            if ~isempty(V0)
                x0 = V0;
%                 assert( norm(obj.eq_lhs*x0(:)-obj.eq_rhs)<1e-6, 'invalid initlization, position constraints not satisfied' );
            else
                x0 = obj.V(:);
                obj.report(2,'Init with %d ARAP iterations  ', obj.initArapIter);
                for arapIter = 1:obj.initArapIter 
                    obj.Tx = obj.T*x0;
                    obj.R = projectRotationMexFast(obj.Tx, obj.dim);
                    x0 = solveConstrainedLS(obj.T, obj.R, obj.eq_lhs, obj.eq_rhs);
                    progBar;
                end
                obj.report(2,'\n');
                x0 = reshape(x0, obj.n_vert, obj.dim);
            end
            % check if solution is orientation preserving
            obj.Tx = obj.T*x0(:);
            [~, ~, obj.flips] = computeFunctionalIsoDistMex(obj.Tx, obj.areas, obj.dim);
            
            
            flipInd2=zeros(obj.n_tri,1);
            pa=obj.Tx;
            J=zeros(obj.n_tri,1);
            for i=1:4:size(pa,1)
                J((i-1)/4+1)=det([pa(i:i+1),pa(i+2:i+3)]);
                if J((i-1)/4+1)<0
                    flipInd2((i-1)/4+1)=1;
                end
            end
            obj.n_flip=sum(flipInd2);
            save dat flipInd2
            k=find(flipInd2);
            if obj.n_flip
                return
            end
            
            
            
            % fix if there are flips
            if obj.flips
                warning('Initialization is not orientation preserving -- projecting on BD')
                % project onto BD
                solver_bd = SolverProjectorBD(obj.F, obj.V, obj.eq_lhs, obj.eq_rhs, obj.proj_bd_K, obj.proj_bd_lb, obj.proj_bd_ub, x0, SolverProjectorModeEnum.Tangent); % setup BD solver
                solver_bd.solve(obj.proj_bd_iter_max, obj.proj_bd_tol_err, obj.proj_bd_verbose); % solve BD projection
                assert(nnz(solver_bd.flips)==0, 'BD projector output has flipped elements');
                x0 = solver_bd.y; % reassign x0
            end
            
            function progBar
                if obj.verbose>=1
                    progressbar(arapIter, obj.initArapIter , 40);
                end
            end
        end
        
        function [varargout] = evaluateFunctional(obj, x, doVal, doGrad, doHess)
            % evaluate
            r=0.3;
            barrier_param=0.000001;
            Vxy=reshape(obj.V,[],1);
            obj.Tx = obj.T*x;
            if doVal||doGrad
                % call mex helps
                [obj.f_val, obj.Tx_grad, obj.flips] = computeFunctionalIsoDistMex(obj.Tx, obj.areas, obj.dim);
                %obj.f_val=obj.f_val+barrier_param*sum(1./(r^2-Vxy.^2));
                if obj.flips
                    warning('negative det');
                    obj.f_val = inf;
                end
                
            end
            % return
            n_arg = 0;
            if doVal
                n_arg = n_arg + 1;
                varargout{n_arg} = obj.f_val;
            end
            if doGrad
                n_arg = n_arg + 1;
                varargout{n_arg} = (obj.Tx_grad'*obj.T)';
                %varargout{n_arg} = (obj.Tx_grad'*obj.T)'+barrier_param*2*Vxy./(r^2-Vxy.^2).^2;
            end
            if doHess
                n_arg = n_arg + 1;
                obj.localHess = computeHessianIsoDistMex(obj.Tx, obj.areas, obj.dim);
                varargout{n_arg} = obj.T'*obj.localHess*obj.T;
            end
        end
        
        function H = setQuadraticProxy(obj)
            wT = spdiag(kron(sqrt(obj.areas), ones(obj.dim^2,1)))*obj.T;
            H = 2*(wT'*wT);
        end
        
        function t_max = getMaxStep(obj, x, p)
            t_max = computeInjectiveStepSizeMex(obj.F,x,p,obj.maxStepTol);
        end
        
        function e = evaluatePerElementEnergy(obj, x)
            % evaluate
            obj.Tx = obj.T*x;
            A = reshape(obj.Tx, obj.dim, obj.dim, []);
            e = nan(obj.n_tri,1);
            for ii = 1:obj.n_tri
                e(ii) = sqrt( sum(sum(A(:,:,ii).^2)) + sum(sum(inv(A(:,:,ii)).^2)) )/sqrt(2*obj.dim);
            end
        end
    end
    
end


