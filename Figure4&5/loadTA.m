function data=loadTA(varargin)
%
% Input options
% name          value           description
% 'NaN'         'interp'        Cubic spline interpolates the NaN points
%                               that appear in the data, truncates the data
%                               on the edges up to the first complete time
%                               trace
%               'zero'          Replaces the NaN values with zeroes
%
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------
%
datafill='interp';

if ~nargin || strncmpi(varargin{1},'gui',3)
    [filename,pathname] = uigetfile('*.csv','Spectrum');%
    if filename==0;
        return; end
    fn=fullfile(pathname,filename);
else
    fn=varargin{1};
end

for i=2:2:nargin-1
    name=varargin{i};
    if ~ischar(name)
        error('Parameter Names Must be Strings.')
    end
    val=varargin{i+1};
    if strncmpi(name,'nan',3), datafill  = val;
    end
end

% throw everything at the file to import it into matlab
try
    rawdata=importdata(fn);
catch %#ok<*CTCH>
    try
        rawdata=load(fn);
    catch
        try
            rawdata=csvread(fn);
        catch
            try
                rawdata = xlsread(fn);
            catch
                try
                    rawdata = dlmread(fn,',');
                catch
                    try
                        rawdata = load(fn,'-ascii');
                    catch
                        try
                            rawdata = readtable(fn,'Delimiter',',','ReadVariableNames',false);
                            rawdata=rawdata{:,:};
                        catch
                            data=[];
                            disp('Trouble loading the data')
                            disp('Check that it is the correct file')
                            disp('Also try truncating the comments off')
                            return;
                        end
                    end
                end
            end
        end
    end
end


%
wavelength = rawdata(2:end,1);
time = rawdata(1,2:end);
rawdata = rawdata(2:end,2:end);

ind=~isnan(wavelength);
wavelength=wavelength(ind);
rawdata=rawdata(ind,:);

% Interpolate the NaN points that might be present in the data
ind=isinf(rawdata);
rawdata(ind)=nan;
if strncmpi(datafill,'i',1)
    % remove the lines that are all NaN from the edges
    ind=isnan(rawdata);
    fi=find(sum(ind,2)~=length(time),1,'first');
    li=find(sum(ind,2)~=length(time),1,'last');
    
    rawdata=rawdata(fi:li,:);
    wavelength=wavelength(fi:li);
    
    % Zero the NaN points on the edges
    ind=isnan(rawdata(1,:));
    rawdata(1,ind)=0;
    ind=isnan(rawdata(end,:));
    rawdata(end,ind)=0;
    
    % Interpolate the NaN values
    dat=rawdata(:);
    dat(isnan(dat)) = interp1(find(~isnan(dat)),...
        dat(~isnan(dat)), find(isnan(dat)), 'spline');
    rawdata=reshape(dat,length(wavelength),length(time));
elseif strncmpi(datafill,'z',1)
    ind=isnan(rawdata);
    rawdata(ind)=0;
end

% round near zero times so when plotted on semilog it doesn't look goofy
% ind=abs(time)<1e-7;
% time(ind)=0;

% output as a structure.
%data = {wavelength(:) time(:) data'};
data.spec=rawdata.';
data.wavelength=wavelength(:);
data.time=time(:);


return