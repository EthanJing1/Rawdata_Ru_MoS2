cd
%% Load Data into the Matlab Environment and perform basic factor analysis
clear %#ok<*ASGLU> 
%#ok<*NASGU>

% loadTA pulls in a dataset, interpolates the NaN points and
% outputs the data into a data structure. 
% If no inputs are provided a file prompt will appear.
% Alternatively the direct file location can be provided such as:
% data = loadTA('C:\Users\Mattk\Desktop\data\ExBox-Per_414nm.csv');

data = loadTA;


% stitchTA will stitch together 2-4 datasets, vis/NIR helios/eos
% it requires an input 'femto' or 'eos' and it will scale the data to the
% provided input
% data = stitchTA('femto');

% This will plot the 2D dataset the function takes up to three inputs.
% The first input needs to be the data cell. The second input designates 
% the time axis scale and can be 'log' or 'linear'. Last is an optional
% parameter which determines the number of contours drawn, default is 30 
% but with really noisy datasets drawing that can be quite slow so less
% contours means it'll draw faster.
figure(1)
plotTA(data,'log',10)

% perform the factor analysis on the data.
figure(2)
EvaluateFactors(data);
%% Truncate Wavelengths
% using logical indices its very easy to truncate the data to produce
% better fits.
% where:
% data.wavelength = wavelength
% data.time = time
% data.spec = spectra, and the first dimension is time and the
% second is wavelength.

% & -and
% | -or
% first define the logical index
ind = data.wavelength>1.75&data.wavelength<1.81;
% then take only the components that would be true in the above index 
data.wavelength = data.wavelength(ind);
data.spec = data.spec(:,ind);

% and replot the results
figure(1)
plotTA(data,'log')

%% Truncate in Time
ind = data.time>7&data.time<7100;
data.time = data.time(ind);
data.spec = data.spec(ind,:);

figure(1)
plotTA(data,'log')

%% First Order Kinetic Model
clear model
% Initial setup of the kinetic model and likely the most important part of
% the fitting. We define buildkineticmodel as a function of the vector t
% where t represents the various unknowns that we are fitting. One can
% supply the initial populations of the different species otherwise the
% algorithm assumes the first species starts at 1 and everything else is 0.

% buildkineticmodel(time,t0,irf,KineticMatrix,population)
model.kinmod=@(t)buildkineticmodel(data,t(1),0.2,...
[-1/t(2) 0; ...
1/t(2) -1/t(3);
],[1 0]);
% Initial guess for values of t

model.initialguess = [0 10 60];
% the labels for t, used for output
model.fitlabels={'Tzero','tau_1','tau_2','tau_3','tau5','tau6','tau7','tau8'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
% model.lb = [-1 0 0 0 ];
%  model.ub = [1 2 200000 50000];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths

% [checkmtime,checkmdecay]= CheckModelSVD(data,model);
selectwavelength=[1.75 1.77 1.79];
CheckModelSW(data,model,selectwavelength);

%% Higher order kinetic model
clear model
clc
% Write out the differential equation where the a's are the population of the
% species and the k's are the rate constants. This equation should be a 
% column vector where each row is a different species.

dadt=@(a,k,t)[-k(1)*a(1);
              k(1)*a(1)-k(2)*a(2);];

% This function builds up the kinetic model, we need to provide the above
% anonymous function dadt. We also need to provide a vector of k's and the 
% initial population.
%
% buildmodeldiff(time,diffeqn,t0,irf,KineticParams,population)
model.kinmod=@(v)buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1 0]);

% The initial guess for the values of v.    
model.initialguess = [0.07 0.3 7.1 79.9];
% also for use later on provide the labels for each parameter.
model.fitlabels={'tau_1','tau_2','t0','irf'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
%model.lb = [0 0 -1 .1];
%model.ub = [10 200  1 .5];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths

CheckModelSVD(data,model)

CheckModelSW(data,model,[ 480 523 650 674 988 1280 ])
%% Decay Associated Model
clear model
% Write out the decay model for the various components. Where each row
% corresponds to a different species. k's are the decay rates or other 
% parameters while t corresponds to time. For an infinite time component
% either provide a component with a fixed long delay or use the following:
% ones(1,length(t))

dt = @(k,t)[exp(-k(1)*t);
          exp(-k(2)*t);
          exp(-k(3)*t)
          
         ];

% With the decay model specified pass it to the function buildmodel which 
% performs the convolution and allows you do determine which parameters
% will vary.
% buildmodel(data,decaymodel,t0,irf,K);
model.kinmod=@(v)buildmodel(data,dt,0,0.6,[1/v(1) 1/v(2) 1/v(3)]);

% initial guess for v
model.initialguess = [1000 10000 100000];
% the labels for v, used for output
model.fitlabels={'irf','tau_2','tau_3','tau_4'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
%model.lb = [-1 0.4 0.9 20 1000];
%model.ub = [0 0.8 2  100 1000000];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths

%CheckModelSVD(data,model)
 
CheckModelSW(data,model,[460 481 519 629 681 1000]);
%% SVD Global Fitting
% The actual fitting and optimization is performed here. The same inputs
% are used as for check kinetic model above, the one exception in the final
% field which is the choice of optimization algorithm. The default is to
% use a Levenberg-Maquardt, 'levmar' algorithm. If that algorithm doesn't
% want to converge a simplex routine can be used, change 'levmar' to 'simplex'
% Simplex will almost always converge to something but is generally much slower.

% SVDfit(data,KineticModel,method)
tic
fit = SVDfit(data,model,'s'); 
toc
% this function will output the results to the command line.
% fitoutput(fit,model);

%% Wavelength Global Fittting
% This routine will fit a single wavelength or range of wavelengths
% simultaneously using the provided model. Input can be discrete values or 
% vectors of values and it will choose the nearest data points.

range = selectwavelength;
tic
fit = SWfit(data,range,model,'l'); 
toc
% this function will output the results to the command line.
% fitoutput(fit,model);
figure;
plot(data.time,fit.Population);legend('species1','species2');
xlabel('time (ps)','FontSize',12,'FontWeight','bold');ylabel('poplulation','FontSize',12,'FontWeight','bold');title('time vs population'); 
%% mesh
negative=fit.Population(:,1)*fit.SpectraFit(1,:);
positive=fit.Population(:,2)*fit.SpectraFit(2,:);
figure
mesh(data.wavelength,data.time,negative)
axis tight
view ([90 -90]);
caxis([-0.02,0.007]);
colorbar;
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 3;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);
ylabel('time (ps)','FontName','Arial','FontSize',12,'FontWeight','bold');
xlabel('probe photon energy(eV)','FontName','Arial','FontSize',12,'FontWeight','bold');
%% 
figure
mesh(data.wavelength,data.time,positive)
axis tight
view ([90 -90]);
caxis([-0.02,0.007]);
colorbar;
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 3;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);
ylabel('time (ps)','FontName','Arial','FontSize',12,'FontWeight','bold');
xlabel('probe photon energy(eV)','FontName','Arial','FontSize',12,'FontWeight','bold');
%% fitting result plotting
figure;

plot(data.time,fit.Population(:,1),'LineWidth',2,'Color','red');
hold on
plot(data.time,fit.Population(:,2),'LineWidth',2,'Color','blue');
axis tight
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 3;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);
legend('Gold','box','off','Location','best');
xlabel('Time (ps)','FontSize',12,'FontWeight','bold','FontName','Arial');
ylabel('Population','FontSize',12,'FontWeight','bold','FontName','Arial');
title('Dynamics','FontSize',14,'FontWeight','bold','FontName','Arial'); 

figure;
plot(data.wavelength,fit.SpectraFit(1,:),'LineWidth',2,'Color','red');
hold on
plot(data.wavelength,fit.SpectraFit(2,:),'LineWidth',2,'Color','blue');
axis tight
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 3;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);
legend('Gold','box','off','Location','best');
xlabel('wavenumber(cm-1)','FontSize',12,'FontWeight','bold','FontName','Arial');
ylabel('Intensity','FontSize',12,'FontWeight','bold','FontName','Arial');
title('Spectra','FontSize',14,'FontWeight','bold','FontName','Arial'); 
sp1=fit.SpectraFit(1,:);
sp2=fit.SpectraFit(2,:);
t=data.time;
w=data.wavelength;
%% plot res and raw
fitdata1=fit.Population(:,1)*fit.SpectraFit(1,:);
fitdata2=fit.Population(:,2)*fit.SpectraFit(2,:);
fitdata=fitdata1+fitdata2;
fitdata1=fitdata*(0.25/0.12);
res=data.spec-fitdata;
figure;mesh(data.time,data.wavelength,fitdata1');
axis tight;
view([0,90]);

colorbar;

ylim([1.69, 1.80]);
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 4;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);
a.YRuler.LineWidth = 4;
xlabel('time (ps)','FontSize',12,'FontWeight','bold','FontName','Arial');ylabel('probe photon energy (eV)','FontSize',12,'FontWeight','bold','FontName','Arial');
title('Ru TMD residual','FontName','Arial','FontSize',14,'FontWeight','bold');
%% 
figure
scatter(fit.time,fit.KineticData(:,1),10,[0.8627, 0.3412, 0.0706],"filled")
hold on
scatter(fit.time,fit.KineticData(:,2),10,[0.9569, 0.8157, 0],"filled")
scatter(fit.time,fit.KineticData(:,3),10,[0.3961, 0.5765, 0.2902],"filled")

plot(fit.time,fit.KineticFit,'LineWidth',2,'color',[0, 0, 0]);
axis tight
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 4;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);
a.YRuler.LineWidth = 4;
xlabel('time (ps)','FontSize',12,'FontWeight','bold','FontName','Arial');ylabel('pump probe intensity','FontSize',12,'FontWeight','bold','FontName','Arial');
legend('1.75eV','1.77eV','1.79eV','box','off', Location='best');
%% shaded area
t11=16.91053;
t12=15.81587;
pppt11=1*exp(-(1/t11)*data.time);
pppt12=1*exp(-(1/t12)*data.time);
pppt1=1*exp(-(1/16.3632)*data.time);

t21=147.4505;
t22=110.4925;
kk=(1/16.3632)/((1/16.3632)-(1/128.9715));
pppt21=kk*exp(-(1/t21)*data.time)-kk*exp(-(1/16.3632)*data.time);
pppt22=kk*exp(-(1/t22)*data.time)-kk*exp(-(1/16.3632)*data.time);
pppt2=kk*exp(-(1/128.9715)*data.time)-kk*exp(-(1/16.3632)*data.time);
figure
plot(data.time,pppt11,data.time,pppt12,'Color','red','LineWidth',0.5)
hold on
plot(data.time,pppt1,'Color','red','LineWidth',2)
plot(data.time,pppt21,data.time,pppt22,'Color','blue','LineWidth',0.5)
plot(data.time,pppt2,'Color','blue','LineWidth',2)
xlim([0 86]);
ylim([0 1]);
a = gca;
set(a,'box','off','color','none');
xticklabels = get(a, 'XTickLabel'); 
set(a, 'XTickLabel', xticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.XRuler.LineWidth = 3;
yticklabels = get(a, 'YTickLabel'); 
set(a, 'YTickLabel', yticklabels, 'FontSize', 10,'FontWeight','bold'); 
a.TickDir = 'out';
a.YRuler.LineWidth = 3;
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',3);
axes(a);

xlabel('Time (ps)','FontSize',12,'FontWeight','bold','FontName','Arial');
ylabel('Population','FontSize',12,'FontWeight','bold','FontName','Arial');
title('Dynamics','FontSize',14,'FontWeight','bold','FontName','Arial'); 

%% Save the output
% Evaluate this and you'll save the species from the fit, the corresponding
% kinetics both experimental and fit and the fit parameters

saveresults(fit)

