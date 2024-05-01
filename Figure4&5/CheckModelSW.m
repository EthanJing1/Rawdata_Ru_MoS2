function [checkMdata,checkMfit] = CheckModelSW(varargin)
% CheckSWModel(data,model,waverange)
%
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------


if nargin<3
    error('Not enough input parameters');
end
data = varargin{1};
spec = data.spec;

model = varargin{2};

kinmodel = model.kinmod;
params = model.initialguess;

decay = kinmodel(params);

% pull out the individual wavelength components
waverange = varargin{3}(:);
wavelength = data.wavelength;
time = data.time;


ind=zeros(length(waverange),1);
vp=0;
for i=1:length(waverange)
    [~,v]=min(abs(wavelength-waverange(i)));
    if vp==v
        v=v+1;
    end
    ind(i)=v;
    vp=v;
end
ind=unique(ind);
%wavel=wavelength(ind);
data=spec(:,ind);

w=(decay'*decay)\(decay'*data);
fit=decay*w;

sf=(decay'*decay)\(decay'*spec);



warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
% plot outputs
figure(1)
clf
if min(time)<0
  ind = time<-min(time);
else
  ind = time<0.01*max(time);
end

s1 = subplot(1,3,1);
plot(time,data,time,fit,'k')
xlim([min(time) max(time(ind))])
ylim([min([data(:);fit(:)]) max([data(:);fit(:)])])
xlabel('Time')

s2 = subplot(1,3,2:3);
plot(time,data,time,fit,'k')
xlim([min(time(~ind)) max(time(~ind))])
ylim([min([data(:);fit(:)]) max([data(:);fit(:)])])
set(gca,'Xscale','log','YTickLabel','')

set(gca, 'ColorOrder',jet(length(waverange)), 'NextPlot', 'replacechildren');
xlabel('Time')
hleg = legend(num2str(waverange,3));
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
checkMdata = data;
checkMfit = fit;
return

