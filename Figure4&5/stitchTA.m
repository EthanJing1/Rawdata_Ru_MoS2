function sdata=stitchTA(varargin)
%
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------
%
% clear
if ~nargin
    scl='eos';
else
    scl=varargin{1};
end
[filename,pathname] = uigetfile({'*.CSV'},'Load spectra','MultiSelect', 'on');%

nc=length(filename);

data=cell(nc,7);
% pull data, sort it
for i=1:nc
    dat=loadTA(fullfile(pathname,filename{i}),'nan','interp');
    data{i,1}=dat.time(:);
    data{i,2}=dat.wavelength(:);
    data{i,3}=dat.spec;
    % collect some info about the spectra
    data{i,4}=min(dat.time);
    data{i,5}=max(dat.time);
    data{i,6}=min(dat.wavelength);
    data{i,7}=max(dat.wavelength);
    if data{i,4}<0 && data{i,4}>-0.1;
        data{i,1}=data{i,1}*1000000;
        data{i,4}=min(data{i,1});
        data{i,5}=max(data{i,1});
    end
end




if nc==2
    data=sortrows(data,5);
    if data{2,5}>data{1,5}*10
        % Stitch time
        sdata=stitchtime(data,scl);
    else
        % stitch wavelength
        sdata=stitchwave(data);
    end    
elseif nc==4
    data=sortrows(data,7);
    sdata=cell(2,7);
    for i=1:2
        dat=stitchtime(data((i:i+1)+(i-1),:),scl);
        sdata{i,1}=dat.time(:);
        sdata{i,2}=dat.wavelength(:);
        sdata{i,3}=dat.spec;
        % collect some info about the spectra
        sdata{i,4}=min(dat.time);
        sdata{i,5}=max(dat.time);
        sdata{i,6}=min(dat.wavelength);
        sdata{i,7}=max(dat.wavelength);
    end
    sdata=stitchwave(sdata);
end

if ~nargout
    data=[0,sdata.time.';sdata.wavelength,sdata.spec.'];
    p=pwd;
    cd(pathname);
    [name,loc] = uiputfile('*.csv','Save Stitched Data As');
    try
        csvwrite(fullfile(loc,name),data)
    catch    
        name=strrep(name, 'csv', 'txt');
        save(fullfile(loc,name),'-ascii', '-double', '-tabs','data');
    end
    cd(p);
end
return

function sdata=stitchtime(data,scl)
        data=sortrows(data,5);
        % first truncate the datasets equally
        wmin = [data{1,6} data{2,6}];
        wmax = [data{1,7} data{2,7}];
        range = [max(wmin) min(wmax)];
        
        ind1 = data{1,2}>=range(1) & data{1,2}<=range(2);
        data{1,2} = data{1,2}(ind1);
        data{1,3} = data{1,3}(:,ind1);
        
        ind2 = data{2,2}>=range(1) & data{2,2}<=range(2);
        data{2,2} = data{2,2}(ind2);
        data{2,3} = data{2,3}(:,ind2);
        % interpolate if the lengths don't match
        if length(data{1,2}) ~= length(data{2,2})
            dim=[length(data{1,2}) length(data{2,2})];
            [~,idx]=sort(dim);
            dat=zeros(length(data{idx(2),1}),dim(idx(1)));
            for j=1:length(data{idx(2),1})
                dat(j,:) = interp1(data{idx(2),2},data{idx(2),3}(j,:),data{idx(1),2},'spline');
            end
            data{idx(2),2}=data{idx(1),2};
            data{idx(2),3}=dat;
        end
        dim = length(data{1,2});
        % determine the overlap and down interpolate the femto to scl
        % the slices
        ind1 = data{1,1}>1000;
        ind2 = data{2,1}>1000 & data{2,1}<data{1,5};
        
        ns=data{2,3}(ind2,:);
        a=zeros(1,dim);
        %dat = zeros(sum(ind2),dim);
        for j=1:dim
            %            dat(:,j) = interp1(data{1,1}(ind1),data{1,3}(ind1,j),data{2,1}(ind2),'spline');
            %            a(j)=(dat(:,j)'*dat(:,j))\(dat(:,j)'*ns(:,j));
            dat = interp1(data{1,1}(ind1),data{1,3}(ind1,j),data{2,1}(ind2),'spline');
            if strncmpi(scl,'eos',3)
                a(j) = (dat'*dat)\(dat'*ns(:,j));
                data{1,3}(:,j) = data{1,3}(:,j)*a(j);
            elseif strncmpi(scl,'femto',3)
                a(j) = (ns(:,j)'*ns(:,j))\(ns(:,j)'*dat);
                data{2,3}(:,j) = data{2,3}(:,j)*a(j);
            end
        end
        % stitch
        ind2 = data{2,1}>data{1,5};
        sdata.time = [data{1,1};data{2,1}(ind2)];
        sdata.wavelength = data{1,2};
        sdata.spec = [data{1,3};data{2,3}(ind2,:)];
        %         plot(data{1,1},data{1,3}(:,100),data{2,1}(ind2),dat(:,100));
return

function sdata=stitchwave(data)
        data=sortrows(data,7);
        % fix the time axis
        tmin = [data{1,4} data{2,4}];
        tmax = [data{1,5} data{2,5}];
        range = [max(tmin) min(tmax)];
        
        ind1 = data{1,1}>=range(1) & data{1,1}<=range(2);
        ind2 = data{2,1}>=range(1) & data{2,1}<=range(2);
        idx=sum(ind1)-sum(ind2);
        if idx>0
            data{2,1} = data{2,1}(ind2);
            data{2,3} = data{2,3}(ind2,:);
            dim = length(data{1,2});
            dat = zeros(sum(ind2),dim);
            for j=1:dim
                dat(:,j) = interp1(data{1,1},data{1,3}(:,j),data{2,1},'spline');
            end
            data{1,1} = data{2,1};
            data{1,3} = dat;
        else
            
            data{1,1} = data{1,1}(ind1);
            data{1,3} = data{1,3}(ind1,:);
            dim = length(data{2,2});
            dat = zeros(sum(ind1),dim);
            for j=1:dim
                dat(:,j) = interp1(data{2,1},data{2,3}(:,j),data{1,1},'spline');
            end
            data{2,1} = data{1,1};
            data{2,3} = dat;
        end
        %plot(data{1,2}(ind1),data{1,3}(60,ind1),data{2,2}(ind2),data{2,3}(60,ind2))
        if data{1,7}<data{2,6}
            % there isn't overlap
            m1 = mean(data{1,3}(:,end-5),2);
            m2 = mean(data{2,3}(:,1:5),2);            
            a=(m2'*m2)\(m2'*m1);
        elseif data{1,7}>data{2,6}
            % there is overlap
            range = [data{2,6}-20 data{1,7}+20];
            ind1 = data{1,2}>=range(1) & data{1,2}<=range(2);
            ind2 = data{2,2}>=range(1) & data{2,2}<=range(2);
            m1 = mean(data{1,3}(:,ind1),2);
            m2 = mean(data{2,3}(:,ind2),2);
            a=(m2'*m2)\(m2'*m1);
        end
        %semilogx(data{1,1},m1,data{1,1},m2)
        sdata.time = data{1,1};
        sdata.wavelength = [data{1,2};data{2,2}];
        sdata.spec = [data{1,3} a*data{2,3}];
        [sdata.wavelength,idx]=sort(sdata.wavelength);
        sdata.spec=sdata.spec(:,idx);
return