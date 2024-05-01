%% Load Data into the Matlab Environment and perform basic factor analysis
clear %#ok<*ASGLU> 

data = loadTA('E:\matlab\tasvd-dev\PdTPP_552nm_fs_VIS.csv');

ind = data.time>-2 &data.time<4;
data.time = data.time(ind);
data.spec = data.spec(ind,:);

ind = data.wavelength>445 & data.wavelength<545 | data.wavelength>560 & data.wavelength<825;
data.wavelength = data.wavelength(ind);
data.spec = data.spec(:,ind);

bkg = repmat(mean(data.spec(data.time<-1,:)),length(data.time),1);

data.spec = - log10(10.^(-data.spec) - 10.^(-bkg) + 1);

figure(1)
pcolor(data.wavelength,data.time, data.spec) ; shading flat; axis tight

%%
%[x,y]=ginput(6);

%f = @(a) sum(( y-a(4)-a(1)*sqrt((a(2).*x.^2-1)./(a(3).*x.^2-1)) ).^2);
%a=fminsearch(f,[1 1 1 1]);
a = [0.4,6,1.98,-0.68]; 
pcolor(data.wavelength,data.time, data.spec) ; shading flat; axis tight
hold all
plot(x,a(1)*sqrt((a(2).*x.^2-1)./(a(3).*x.^2-1)+a(4)))
hold off
%%
i = 1;
plot(data.time, data.spec(:,i))


%% First Order Kinetic Model
clear model
% Initial setup of the kinetic model and likely the most important part of
% the fitting. We define buildkineticmodel as a function of the vector t
% where t represents the various unknowns that we are fitting. One can
% supply the initial populations of the different species otherwise the
% algorithm assumes the first species starts at 1 and everything else is 0.

% buildkineticmodel(time,t0,irf,KineticMatrix,population)
model.kinmod=@(t)buildkineticmodel(data,t(1),t(2),...
    [-1/t(3) 0;
     1/t(3) -1/t(4);],[1 0]);

% Initial guess for values of t
model.initialguess = [0 0.3  7.1 79.9];
% the labels for t, used for output
model.fitlabels={'t0','irf','tau_1','tau_2'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
%model.lb = [0 0 -1 .1];
%model.ub = [10 200  1 .5];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths
tic
CheckModelSVD(data,model)
toc
%CheckModelSW(data,model,[529 704 951])
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

% CheckModelSW(data,model,[529 704 951])
%% Decay Associated Model
clear model
% Write out the decay model for the various components. Where each row
% corresponds to a different species. k's are the decay rates or other 
% parameters while t corresponds to time. For an infinite time component
% either provide a component with a fixed long delay or use the following:
% ones(1,length(t))

dt=@(k,t)[exp(-k(1)*t);
          exp(-k(2)*t)];

% With the decay model specified pass it to the function buildmodel which 
% performs the convolution and allows you do determine which parameters
% will vary.
% buildmodel(data,decaymodel,t0,irf,K);
model.kinmod=@(v)buildmodel(data,dt,v(1),v(2),[ 1/v(3) 1/v(4)]);

% initial guess for v
model.initialguess = [0.12 0.3 6 60];
% the labels for v, used for output
model.fitlabels={'t0','irf','tau_2','tau_3'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
%model.lb = [0 0 -1 .1];
%model.ub = [10 200  1 .5];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths

CheckModelSVD(data,model)
 
%CheckModelSW(data,model,[529 704 951])
%% SVD Global Fitting
% The actual fitting and optimization is performed here. The same inputs
% are used as for check kinetic model above, the one exception in the final
% field which is the choice of optimization algorithm. The default is to
% use a Levenberg-Maquardt, 'levmar' algorithm. If that algorithm doesn't
% want to converge a simplex routine can be used, change 'levmar' to 'simplex'
% Simplex will almost always converge to something but is generally much slower.

% SVDfit(data,KineticModel,method)
tic
fit = SVDfit(data,model,'l'); 
toc
% this function will output the results to the command line.
% fitoutput(fit,model);

%% Wavelength Global Fittting
% This routine will fit a single wavelength or range of wavelengths
% simultaneously using the provided model. Input can be discrete values or 
% vectors of values and it will choose the nearest data points.

range = [600 700];
tic
fit = SWfit(data,range,model,'l'); 
toc
% this function will output the results to the command line.
% fitoutput(fit,model);

%% Save the output
% Evaluate this and you'll save the species from the fit, the corresponding
% kinetics both experimental and fit and the fit parameters

saveresults(fit)

