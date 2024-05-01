function fit = fitoutput(fit,model)
% fitoutput(labels,fitparameters,statistics)

time = fit.time;
kfit = fit.KineticFit;
kin = fit.KineticData;

wavelength = fit.wavelength;
sf = fit.SpectraFit;

labels = model.fitlabels;
params = fit.FitParams;
std = fit.Std;
var = fit.MeanSquareError;

warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
% plot outputs
% figure(1)
% clf
% subplot(2,1,1)
% plot(time,kfit,time,kin)
% xlim([0 max(time)])
% set(gca,'Xscale','log')
%
% subplot(2,1,2)
% plot(time,kfit,time,kin)
% xlim([min(time) max(time)])
% drawnow
% xlabel('Time')

figure(1)
clf
if min(time)<0
  ind = time<-min(time);
else
  ind = time<0.01*max(time);
end

s1 = subplot(1,3,1);
plot(time,kin,time,kfit,'k')
xlim([min(time) max(time(ind))])
ylim([min([kfit(:);kin(:)]) max([kfit(:);kin(:)])])
xlabel('Time')

s2 = subplot(1,3,2:3);
plot(time,kin,time,kfit,'k')
xlim([min(time(~ind)) max(time(~ind))])
ylim([min([kfit(:);kin(:)]) max([kfit(:);kin(:)])])
set(gca,'Xscale','log','YTickLabel','')


xlabel('Time')
if isfield(fit,'waverange')
  hleg = legend(num2str(round(fit.waverange,2)));
  set(hleg,'Location','best')
else
  l = size(kfit,2);
  leg = 'Species ';
  leg = repmat(leg,l,1);
  leg = [leg,num2str((1:l)')];
  hleg = legend(leg);
  set(hleg,'Location','best')
end

p1 = get(s1,'position');
p2 = get(s2,'position');
set(s2,'position',[p1(1)+p1(3) p1(2) p2(3)+(p2(1)-(p1(1)+p1(3))) p1(4)])
linkaxes([s1 s2],'y')
drawnow

figure(2)
clf
plot(wavelength,sf,wavelength,zeros(length(wavelength),1))
xlim([min(wavelength) max(wavelength)])
xlabel('Wavelength (nm)')

l = size(sf,1);
leg = 'Species ';
leg = repmat(leg,l,1);
leg = [leg,num2str((1:l)')];
hleg = legend(leg);
set(hleg,'Location','best')

drawnow

nout = min([length(labels) length(params)]);
str=['\n',datestr(now),'\n','MSE       =  ',num2str(var,'%6.4e'),'\n'];
disp(' ')
disp(datestr(now))
disp(' ')
disp(['MSE       =  ',num2str(var,'%6.4e')])
for i = 1:nout
  str = [str, labels{i} '  =  ',num2str(params(i),'%6.4f'),' +- ',num2str(std(i),'%6.4e'),'\n']; %#ok<*AGROW>
  disp([labels{i} '  =  ',num2str(params(i),'%6.4f'),' +- ',num2str(std(i),'%6.4e')])
end
str = [str,'\n',func2str(model.kinmod),'\n'];
disp(' ')
fit.str = str;
return
% if we want to get fancy we can use this
% text('Interpreter','latex',...
%  'String','$$\int_0^x\!\int_y dF(u,v)$$')