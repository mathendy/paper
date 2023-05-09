function triEn = energyForEverySingleTriangle(fz,gz,areas)
%ENERGYFOREVERYSINGLETRIANGLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
fz2=real(fz).*real(fz)+imag(fz).*imag(fz);
gz2=real(gz).*real(gz)+imag(gz).*imag(gz);
triEn=(fz2+gz2).*(1+1./(gz2-fz2).^2).*areas;
end

