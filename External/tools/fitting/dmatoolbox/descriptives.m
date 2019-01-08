function varargout = descriptives(fulldata,fn)
%DESCRIPTIVES  Compute descriptive statistics for two-choice RT data
%   TABLE = DESCRIPTIVES(DATA), where DATA is a valid N-by-3 or N-by-2 data
%   set, returns a cell matrix with some descriptive statistics of the RT
%   and choice data.
%
%   DESCRIPTIVES(DATA) prints the table to the screen.
%
%   DESCRIPTIVES(DATA,FILENAME), prints the table to the file FILENAME.
%
%   See also QPPLOT and EDFCDF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Fixed typo that made function crash if input was N-by-2.

if nargin<1
    error('DMAT:descriptives:notEnoughInputs',...
        'DESCRIPTIVES requires at least one input variable.')
end

if nargin>1
    fid = fopen(fn,'w');
else
    fid = 1;
end

if isvaliddataset(fulldata,3)
    conds = unique(fulldata(:,1));
    ncond = length(conds);
elseif isvaliddataset(fulldata,2)
    ncond = 1;
    conds = 1;
    fulldata = [ones(size(fulldata,1),1) fulldata];
else
    error('DMAT:descriptives:invalidDataSet',...
        'This is not a valid data set.')
end

dc = {};

for ctr=0:ncond
    if ctr
        data = fulldata(fulldata(:,1)==conds(ctr),:);
        conname=sprintf('Condition %i',conds(ctr));
    else
        data = fulldata;
        conname='Across conditions';
    end

    N = length(data(:,2));
    nx1 = sum(data(:,2));   nx0 = N-nx1;
    px1 = nx1/N;            px0 = 1-px1;
    stdp = sqrt(px1*px0/N);

    nx1 = num2str(nx1,'%8.0f');
    px1 = num2str(px1,'%8.5f');
    nx0 = num2str(nx0,'%8.0f');
    px0 = num2str(px0,'%8.5f');
    stdp = num2str(stdp,'%8.5f');
    N = num2str(N,'%8.0f');

    RT = data(:,3);
    RT0=RT(data(:,2)==0);
    RT1=RT(data(:,2)==1);

    warning('off','MATLAB:divideByZero')
    avgRT  = num2str(mean(RT),'%8.5f');   stdRT  = num2str(std(RT),'%8.5f');
    maxRT  = num2str(max(RT),'%8.5f');    minRT  = num2str(min(RT),'%8.5f');
    avgRT0 = num2str(mean(RT0),'%8.5f');  avgRT1 = num2str(mean(RT1),'%8.5f');
    maxRT0 = num2str(max(RT0),'%8.5f');   maxRT1 = num2str(max(RT1),'%8.5f');
    minRT0 = num2str(min(RT0),'%8.5f');   minRT1 = num2str(min(RT1),'%8.5f');
    stdRT0 = num2str(std(RT0),'%8.5f');   stdRT1 = num2str(std(RT1),'%8.5f');
    warning('on','MATLAB:divideByZero')

    co = {conname       ,''     ,''              ,''
        'Accuracy'      ,''     ,''              ,''
        'P(X=0)'        ,px0    ,'P(X=1)'        ,px1
        'N(X=0)'        ,nx0    ,'N(X=1)'        ,nx1
        'STD(P[X=0])'   ,stdp   ,'STD(P[X=1])'   ,stdp
        'Total N'       ,N      ,''              ,''
        ''              ,''     ,''              ,''
        'Response time' ,''     ,''              ,''
        'Mean(RT)'      ,avgRT  ,'STD(RT)'       ,stdRT
        'Max(RT)'       ,maxRT  ,'Min(RT)'       ,minRT
        'Mean(RT[X=0])' ,avgRT0 ,'STD(RT[X=0])'  ,stdRT0
        'Max(RT[X=0])'  ,maxRT0 ,'Min(RT[X=0])'  ,minRT0
        'Mean(RT[X=1])' ,avgRT1 ,'STD(RT[X=1])'  ,stdRT1
        'Max(RT[X=1])'  ,maxRT1 ,'Min(RT[X=1])'  ,minRT1
        ''              ,''     ,''              ,''
        ''              ,''     ,''              ,''    };
    dc = [dc;co];
end

if nargout
    varargout{1} = dc;
else        
    for row = 1:size(dc,1)
        fprintf(fid,'%15s',dc{row,1});
        fprintf(fid,'%10s',dc{row,2});
        fprintf(fid,'%20s',dc{row,3});
        fprintf(fid,'%10s\n',dc{row,4});
    end
    if fid~=1
        fclose(fid);
    end
end