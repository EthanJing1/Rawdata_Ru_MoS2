function plotTA(varargin)
clf
data=varargin{1};
scale=varargin{2};

spec = data.spec;
wavelength = data.wavelength(:).';
time = data.time(:);

warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
subplot(2,1,1)

if strncmpi(scale,'log',3)
    ind=time>0;
    pcolor(wavelength,time(ind),spec(ind,:)); shading flat; 
    %contourf(wavelength,time(ind),spec(ind,:),cont,'LineStyle','none')
    set(gca,'Yscale','log')
    
elseif strncmpi(scale,'lin',3)
    pcolor(wavelength,time,spec); shading flat; 
    %contourf(wavelength,time,spec,cont,'LineStyle','none')
end
ylabel('time')
xlim([min(wavelength) max(wavelength)])

subplot(2,1,2)
in=time>0;

ind=find(in,1)+0:ceil(length(time(in))/6):length(time); 

ind=[3 ind];
set(gca, 'ColorOrder',jet(length(ind)), 'NextPlot', 'replacechildren');
plot(wavelength,spec(ind,:))

legend(num2str(time(ind)))
set(gca,'yticklabel','')
xlabel('Wavelength')
xlim([min(wavelength) max(wavelength)])

return
