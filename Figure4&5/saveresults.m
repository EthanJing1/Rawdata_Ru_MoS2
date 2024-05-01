function saveresults(fit) %#ok<*INUSL,*INUSD>
% outputs the
[name,loc] = uiputfile({'*.txt';'*.csv'});

wavelength = fit.wavelength;
sfit = fit.SpectraFit;
sfit = [wavelength(:), sfit.']; %#ok<*NASGU>

fitpar = fit.str;

if isfield(fit,'waverange')
  
  time = fit.time;
  range = fit.waverange;
  
  kfit = fit.KineticFit;
  kfit = [0,range(:)';time(:),kfit];
  
  kin = fit.KineticData;
  kdat = [0,range(:)';time(:) kin];
  
else
  time = fit.time;
  kfit = fit.KineticFit;
  kfit = [time(:),kfit];
  
  kin = fit.KineticData;
  kdat = [time(:) kin];
  
end


[~,name,ext]=fileparts(name);
switch ext
  case '.txt'
    % save the data in tab delimitted
    save(fullfile(loc,[name,'_spectra.txt']),'-ascii', '-double', '-tabs','sfit');
    save(fullfile(loc,[name,'_kinfit.txt']),'-ascii', '-double', '-tabs','kfit');
    save(fullfile(loc,[name,'_kindat.txt']),'-ascii', '-double', '-tabs','kdat');
    
    fid = fopen(fullfile(loc,[name,'_fitpar.txt']),'wt');
    %fwrite(fid, fitpar, 'uchar');
    fprintf(fid, fitpar);
    fclose(fid);
  case '.csv'
    csvwrite(fullfile(loc,[name,'_spectra.csv']),sfit);
    csvwrite(fullfile(loc,[name,'_kinfit.csv']),kfit);
    csvwrite(fullfile(loc,[name,'_kindat.csv']),kdat);
    
    fid = fopen(fullfile(loc,[name,'_fitpar.txt']),'wt');
    %fwrite(fid, fitpar, 'uchar');
    fprintf(fid, fitpar);
    fclose(fid);
end
return