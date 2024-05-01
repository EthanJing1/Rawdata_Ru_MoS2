function [spec1,spec2]=matchTA(varargin)
% [spec1,spec2]=matchTA(spec1,spec2,match)
% match - 'wavelength'
%         'time'
%         'both' (default if no input)
%
% Matches the spectral and time ranges of two TA spectra for ease of
% subtraction
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------
spec1=varargin{1};
spec2=varargin{2};
if nargin>2
  match=varargin{3};
else
  match='both';
end

if strcmpi(match,'wavelength') || strcmpi(match,'both')
  wlmin = max([min(spec1.wavelength) min(spec2.wavelength)]);
  wlmax = min([max(spec1.wavelength) max(spec2.wavelength)]);
  
  ind = spec1.wavelength>=wlmin & spec1.wavelength<=wlmax;
  spec1.wavelength=spec1.wavelength(ind);
  spec1.spec=spec1.spec(:,ind);
  
  ind = spec2.wavelength>=wlmin & spec2.wavelength<=wlmax;
  spec2.wavelength=spec2.wavelength(ind);
  spec2.spec=spec2.spec(:,ind);
  
  if length(spec1.wavelength)~=length(spec2.wavelength)
    [~,I]=min([length(spec1.wavelength) length(spec2.wavelength)]);
    if I==1
      dat=zeros(length(spec2.time),length(spec1.wavelength));
      for i=1:length(spec2.time)
        dat(i,:)=interp1(spec2.spec(i,:),spec2.wavelength(:)',spec1.wavelength(:)','spline');
      end
      spec2.spec=dat;
      spec2.wavelength=spec1.wavelength;
    else
      dat=zeros(length(spec1.time),length(spec2.wavelength));
      for i=1:length(spec1.time)
        dat(i,:)=interp1(spec1.wavelength(:)',spec1.spec(i,:),spec2.wavelength(:)','spline');
      end
      spec1.spec=dat;
      spec1.wavelength=spec2.wavelength;
    end
  end
end

if strcmpi(match,'time') || strcmpi(match,'both')
  tmin = max([min(spec1.time) min(spec2.time)]);
  tmax = min([max(spec1.time) max(spec2.time)]);
  
  ind = spec1.time>=tmin & spec1.time<=tmax;
  spec1.time=spec1.time(ind);
  spec1.spec=spec1.spec(ind,:);
  
  ind = spec2.time>=tmin & spec2.time<=tmax;
  spec2.time=spec2.time(ind);
  spec2.spec=spec2.spec(ind,:);
  
  if length(spec1.time)~=length(spec2.time)
    [~,I]=min([length(spec1.time) length(spec2.time)]);
    if I==1
      dat=zeros(length(spec1.time),length(spec2.wavelength));
      for i=1:length(spec2.wavelength)
        dat(:,i)=interp1(spec2.time(:)',spec2.spec(:,i),spec1.time(:)','spline');
      end
      spec2.spec=dat;
      spec2.time=spec1.time;
    else
      dat=zeros(length(spec2.time),length(spec1.wavelength));
      for i=1:length(spec1.wavelength)
        dat(:,i)=interp1(spec1.time(:)',spec1.spec(:,i),spec2.time(:)','spline');
      end
      spec1.spec=dat;
      spec1.time=spec2.time;
    end
  end
end

return
