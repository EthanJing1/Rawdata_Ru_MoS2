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
%data.time=data.time*1000 %convert to ns for nsTA

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
ind = data.wavelength>431&data.wavelength<808 | data.wavelength>870 &data.wavelength<1800;
% then take only the components that would be true in the above index 
data.wavelength = data.wavelength(ind);
data.spec = data.spec(:,ind);

% and replot the results
figure(1)
plotTA(data,'log')

%% Truncate in Time
ind = data.time>-1.7 & data.time <7000; 
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
model.kinmod=@(t)buildkineticmodel(data,t(1),0.4,...
[-1/t(2),0,0;
1/t(2),-1/t(3),0;
0,0,0],[0.75,0,0])+buildkineticmodel(data,t(1),0.4,...
[-1/t(2),0,0;
1/t(2),-1/t(4),0;
0,1/t(4),-1/t(5)],[0.25,0,0]);
% Initial guess for values of t
model.initialguess = [0 55 500 10 10000];
% the labels for t, used for output
model.fitlabels={'tau_1','tau_2','tau3','tau4','tau5'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
%model.lb = [0 ];
%model.ub = [10 ];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths

%[checkMtime,checkMdecay] = CheckModelSVD(data,model);

[checkMdata,checkMfit]=CheckModelSW(data,model,[470 570 634 725 848 931]);
%RMSE = sqrt(sum(sum((checkMdata-checkMfit),2).^2)/267);

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

dt = @(k,t)[exp(-k(1)*t);
          exp(-k(2)*t);
          exp(-k(3)*t)];

% With the decay model specified pass it to the function buildmodel which 
% performs the convolution and allows you do determine which parameters
% will vary.
% buildmodel(data,decaymodel,t0,irf,K);
model.kinmod=@(v)buildmodel(data,dt,v(1),v(2),[ 1/v(3) 1/v(4) 1/v(5)]);

% initial guess for v
model.initialguess = [0 0.6 0.6 11 21000];
% the labels for v, used for output
model.fitlabels={'t0','irf','tau_2','tau_3'};

% To constrain the fitting provide lower and upper bounds for all the 
% parameters, this input in optional.  
%model.lb = [0 0 -1 .1];
%model.ub = [10 200  1 .5];

% check the model and adjust the initial guessusing either SVD or applying 
% it at single/multiple wavelengths

%[a,b] = CheckModelSVD(data,model)
 
CheckModelSW(data,model,[450 570 607 634 660 725])
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

range = [470 570 634 725 848 931];
tic
fit = SWfit(data,range,model,'s'); 
toc
% this function will output the results to the command line.
% fitoutput(fit,model);

%% Save the output
% Evaluate this and you'll save the species from the fit, the corresponding
% kinetics both experimental and fit and the fit parameters

saveresults(fit)

