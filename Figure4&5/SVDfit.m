function fit = SVDfit(varargin)
% bestfit=SVDfit(data,model,algorithm)
%
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------


% check and assign the inputs found in varargin
if nargin<2
    disp('Not enough input parameters');return
end
data = varargin{1};
fit = data;
data=data.spec;

model = varargin{2};

% the number of principal components to use is hardcoded to be the number
% of chemical spectral species present from the kinetics
decay=model.kinmod(model.initialguess);
numsp=size(decay,2);

model = fixmodel(model);
guess=model.initialguess;

if nargin==3
    algorithm=varargin{3};
else
    algorithm='levmar';
end

% Perform SVD on absorption data
[U,S,V] = svd(data);

% Truncate the SVD output based on the number of components we wish to use
V = V(:,1:numsp);
S = S(1:numsp,1:numsp);
U = U(:,1:numsp);

% Just to make sure things scale appropriately later, we'll fold S into V
sv=S*V';

% res = data-U*S*V';

% Fit the kinetic model using a linear combination of the vectors in U by
% exploiting the orthogonality of the columns in U.
if strncmp(algorithm,'l',1)
    % Levenberg-Maquardt optimization Converges very fast, but can be
    % sensitive to the starting guess
    res=@(params)fitkin(params,model,U);
    [bestfit, ~, iter] = LMFnlsq(res,guess,'MaxIter', 1000,'FunTol',1e-20,'XTol',1e-10);
    %    bestfit  = smarquardt(res, guess);
    if iter==1 || iter ==-100
        disp('Bad initial guess or poor model');return
    end
elseif strncmp(algorithm,'s',1)
    % Built in Matlab simplex routine can be very slow but will almost
    % always find a minimum.
    
    res=@(params)fitkinSSE(params,model,U);
    
    options=optimset('MaxFunEvals',3000, 'MaxIter', 3000,'TolFun',1e-10,'TolX',1e-10);
    [bestfit,~,flag] = fminsearch(res,guess,options);
    if flag<1
        disp('Bad initial guess or poor model');return
    end
end

bestfit = model.normfun(bestfit);
% sort of a work around at the moment for dealing with the constraints
if isfield(model,'lb') && isfield(model,'ub')
    model = rmfield(model,['lb'; 'ub']);
end
res=@(params)fitkin(params,model,U);
% Calculate the fit based on minimum obtained through fitting.
[~,R,kfit]=res(bestfit);
[mse,std]=lmstats(res,bestfit);

kin=U*R;
spec=R\sv;
sf=(kfit'*kfit)\(kfit'*data);

fd=kfit*sf;
residual=data-fd;
% store the results in the fit structure
fit.FitParams = bestfit;
fit.KineticFit = kfit;
fit.KineticData = kin;
fit.SpectraFit = sf;
fit.Residual = residual;
fit.MeanSquareError = mse;
fit.Std = std;

fit = fitoutput(fit,model);

return

function model = fixmodel(model)
if isfield(model,'lb') && isfield(model,'ub')
    lb = model.lb;
    ub = model.ub;
    guess = model.initialguess;
    % rescale the guess to take into account the bounds
    x = [lb;ub];
    range = max(x) - min(x);
    mx = mean(x);
    model.initialguess = 2*(guess - mx )./range;
    model.normfun = @(y) range(:).*y(:)/2+mx(:);
else
    model.normfun = @(y) y(:);
end
return


function [res,r,decay]=fitkin(params,model,U)
if isfield(model,'lb') && isfield(model,'ub')
    ind = abs(params)>1;
    penalty = sum(params(ind).^2);
    params = model.normfun(params);
else
    penalty = 0;
end

kinmodel = model.kinmod;
decay = kinmodel(params);

% Linear combination of U used to fit the decay
r=U'*decay;
res=U-decay/r;

% For Levenberg-Maquardt optimization we need a column of the residuals to
% to calculate the jacobian
res=res(:)+penalty;

return

function [res,r,decay]=fitkinSSE(params,model,U)

if isfield(model,'lb') && isfield(model,'ub')
    ind = abs(params)>1;
    penalty = sum(params(ind).^2);
    params = model.normfun(params);
else
    penalty = 0;
end

kinmodel = model.kinmod;
decay = kinmodel(params);

% Linear combination of U used to fit the decay

r=U'*decay;
res=U-decay/r;

% Sum of the Squares of residual
% for simplex
res=res(:)'*res(:)+penalty;

return
