

%%
clear
data = loadMFE;
%% Baseline

ind = data.time<-300;
for i = 1:length(data.wavelength)
  bl = mean(data.spec(ind,i));
  data.spec(:,i) = data.spec(:,i)-bl;
end


%% Determine the zerofield points
ind = data.wavelength <= 5;
index = 1:length(data.wavelength);

% save the index later to account for degradation
data.zfind = index(ind);
data.index = index(~ind);
data.field = data.wavelength(:,~ind);
data.wavelength = index;

%% explore the data
i = 47;
plot(1:length(data.time),data.spec(:,i));


%% Decay Associated Model
clear model

dt=@(k,t)[exp(-abs(k(1))*t);
          exp(-abs(k(2))*t);
          ones(1,length(t));];

model.kinmod=@(v)buildmodel(data,dt,v(1),v(2),[ 1/2 1/v(3) ]);

% initial guess for v
model.initialguess = [-246 3.9 11 ];
% the labels for v, used for output
model.fitlabels={'t0','tau_1','tau_2'};



CheckModelSW(data,model,1)
%% fit the data individually, correct for decomposition and plot MFE

fits = zeros(length(data.wavelength),length(model.initialguess));
ncomps = size(model.kinmod(model.initialguess),2);
amp = zeros(ncomps,length(data.wavelength));

for i = 1:length(data.wavelength)
  range = data.wavelength(i);
  
  fit = SWfit(data,range,model,'l');
  fits(i,:) = fit.FitParams;
  amp(:,i) = fit.SpectraFit(:,i);
end

% correct for decomposition
bl = zeros(size(amp));
for i = 1:ncomps
  [p,S,mu] = polyfit(data.zfind,amp(i,data.zfind),1);
  bl(i,:) = polyval(p,data.wavelength,S,mu);
end
MFE = amp./bl;
figure(1)
plot(data.field,MFE(2:3,data.index))

figure(2)
plot(data.field,fits(data.index,3))

