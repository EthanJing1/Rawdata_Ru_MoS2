function [CheckMtime,CheckMdecay]= CheckModelSVD(data,model)

%
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------



% Pull the individual components out of the data cell


wavelength = data.wavelength;
time = data.time;
data = data.spec;

kinmodel = model.kinmod;
params=model.initialguess;
% We need to provide an initial guess for the values of t.    


decay=kinmodel(params);
numsp=size(decay,2);

% Perform SVD on absorption data
[U,S,V] = svd(data);

% Truncate the SVD output based on the number of components we wish to use
V = V(:,1:numsp);
S = S(1:numsp,1:numsp);
U = U(:,1:numsp);

% Just to make sure things scale appropriately later we'll fold S into V
sv=S*V';

% calculate a rotation matrix such that U will optimally fit the kinetic 
% trace. We're exploiting the orthogonality of U to calculate R.
R=U'*decay;
kin=U*R;
res=kin-decay;

spec=R\sv;
sf=decay\data;

warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
% plot outputs
% figure(1)
% clf
% subplot(2,1,1)
% plot(time,decay,time,kin)
% xlim([0 max(time)])
% set(gca,'Xscale','log')
% 
% subplot(2,1,2)
% plot(time,decay,time,kin)
% xlim([min(time) -min(time)])
% drawnow
% xlabel('Time')
%  
figure(1)
clf
if min(time)<0
  ind = time<-min(time);
else
  ind = time<0.01*max(time);
end

s1 = subplot(1,3,1);
plot(time,decay,time,kin)
xlim([min(time) max(time(ind))])
ylim([min([decay(:);kin(:)]) max([decay(:);kin(:)])])
xlabel('Time')

s2 = subplot(1,3,2:3);
plot(time,decay,time,kin)
xlim([min(time(~ind)) max(time(~ind))])
ylim([min([decay(:);kin(:)]) max([decay(:);kin(:)])])
set(gca,'Xscale','log','YTickLabel','')

xlabel('Time')

l = size(decay,2);
leg = 'Species ';
leg = repmat(leg,l,1);
leg = [leg,num2str((1:l)')];
hleg = legend(leg);
set(hleg,'Location','best')

p1 = get(s1,'position');
p2 = get(s2,'position');
set(s2,'position',[p1(1)+p1(3) p1(2) p2(3)+(p2(1)-(p1(1)+p1(3))) p1(4)])
linkaxes([s1 s2],'y')
drawnow

figure(2)
clf
plot(wavelength,sf)
xlim([min(wavelength) max(wavelength)])
xlabel('wavelength')

l = size(sf,1);
leg = 'Species ';
leg = repmat(leg,l,1);
leg = [leg,num2str((1:l)')];
hleg = legend(leg);
set(hleg,'Location','best')

drawnow
CheckMtime = time;
CheckMdecay = decay;
return

