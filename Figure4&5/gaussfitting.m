clear 
data = loadTA;

plotTA(data,'linear')
%% pull out stokes or antistokes

ind = data.wavelength >-300 & data.wavelength<-90;

s.wavelength = data.wavelength(ind);

s.spec = data.spec(:,ind);
s.time = data.time;

plotTA(s,'linear')

%%
i = 1;
[~,fit] = fitgauss([0.0336 -215 22 0.0066 -231 22],s.spec(i,:),s.wavelength,5);

plot(s.wavelength,s.spec(i,:),s.wavelength,fit)

%%
i = 3;
res = @(a)fitgauss([a(1) a(2) a(3) a(4) a(5) a(6)],s.spec(i,:),s.wavelength,5);
A = [0.0336 -215 22 0.0066 -231 22];
par = fminsearch(res,A);

[~,fit,g,bl] = fitgauss(par,s.spec(i,:),s.wavelength,5);

plot(s.wavelength,s.spec(i,:),s.wavelength,fit)
%plot(s.wavelength,s.spec(i,:)-bl,s.wavelength,g)
%%


A = [0.188,-218.,21.,-0.00187,-321.74,4.079];

par = zeros(6,length(s.time));
for i = 1:length(s.time)
  res = @(a)fitgauss([a(1) a(2) a(3) a(4) a(5) a(6)],s.spec(i,:),s.wavelength,5);
  par(:,i) = fminsearch(res,A);
end










