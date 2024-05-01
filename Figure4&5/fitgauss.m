function [ssr,fit,g,bl] = fitgauss(var,data,wavelength,p)

n = length(var);
n = n/3;
g = zeros(n,length(wavelength));
for i=1:n
  g(i,:) = var(1+(i-1)*3)*gaussian(wavelength,var(2+(i-1)*3),abs(var(3+(i-1)*3)) );
end

res = data - sum(g,1);

[p,S,mu] = polyfit(wavelength',res,p);
bl = polyval(p,wavelength',S,mu);

res = res-bl;
fit = sum(g,1) + bl;
ssr = sum(res.^2);