function data = loadMFE
% LOAD MFE FILES
% Each kinetic is in a separate file.  Load all of the files corresponding to the selected dataset

%open load file dialog
[filename,pathnameload] = uigetfile('*.*','Select one kinetic to load dataset, or select all kinetics to process','MultiSelect', 'on');
if pathnameload==0;
  return; end
if ~iscell(filename) %with multiselect, filename is a cell, if it's not either user cancelled load or only picked one file
  
  %user only picked one file: get file base name and load all files matching pattern
  basename = filename(1:end-3); %truncate ### from end of filename
  clear filename;
  for i=1:999
    if i > 99
      zstr = ''; %ie. file is 'basename100'
    elseif i > 9
      zstr = '0'; %file is 'basename010'
    else
      zstr = '00'; %file is 'basename001'
    end;
    
    testname = strcat(basename, zstr, int2str(i));
    
    if exist(strcat(pathnameload, testname), 'file') == 2 %if the name exists as a file add it to the filename array
      filename{i} = testname;
    else
      clear testname zstr;
      break;
    end;
  end;
end;

%load the files one at a time, extract field and other parameters, and build arrays of data
for i = 1:size(filename, 2)
  fid = fopen(strcat(pathnameload, filename{i}), 'r'); %open file, get ID
  rawdata = textscan(fid,'%s%s');
  
  
  
  
  if i == 1 %if first file, populate time axis of data array
    % determine length of data
    indx = findstring(rawdata{1},'[end]');    
    data.spec = zeros(size(filename, 2), size(rawdata{1, 1},1) - indx); %get size of time axis plus one and make data array
    data.time = str2double(rawdata{1, 1}(indx+1:end));
    data.wavelength = zeros(size(filename, 2),1);
  end;
  
  ind = findstring(rawdata{1},'Current');
  %put magnetic field into its axis - will sort later
  dstring = strsplit(rawdata{1, 2}{ind},':'); %split out the current magnetic field from the line containing this data
  data.wavelength(i) = str2double(dstring{2}(1:end-1));%trim off the units 'G' from end and convert to number
  
  %fill in kinetic data
  data.spec(i,:) = str2double(rawdata{1, 2}(indx+1:end));
  
  %close out file
  fclose(fid);
end;

data.spec = data.spec';
data.wavelength = data.wavelength(:).';
data.time = data.time(:);

fclose all; %close all files, probably not needed
return

function indx = findstring(tab,str)
for indx=1:length(tab)
    matches = strfind(tab{indx},str);
    if (~isempty(matches))
        return;
    end;
end;
indx=0;
return;