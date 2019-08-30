function [Proxh, Proxv] = prox_L12(Ph,Pv,gamma)

Proxh = zeros(size(Ph));
Proxv = zeros(size(Pv));

temp=(Ph.^2+Pv.^2).^(1/2);
ind=find(temp>gamma);
Proxh(ind)=(1-gamma./temp(ind)).*Ph(ind);
Proxv(ind)=(1-gamma./temp(ind)).*Pv(ind);

end
