function outputtoxls(data,output,fn)
%OUTPUTTOXLS  Export DMAT output to structured Excel spreadsheet
%   OUTPUTTOXLS(DATA,OUTPUT,FN), where DATA is an N-by-3 data matrix,
%   OUTPUT is a DMAT output structure (or an array of such), and FN is the
%   filename of an Excel spreadsheet.
%
%   If you call this function from the base workspace on a PC, the Excel
%   file will automatically be opened.
%
%   WARNING: Existing files by the name of FN will be overwritten without
%   further warning.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Added catch if Excel not installed on PC.
%  Edit 0.4: Fixed typo in exporting modelfittable.

%% Input checking
if nargin<3
    error('DMAT:outputtoxls:notEnoughInputs',...
        'OUTPUTTOXLS requires three input variables.')
end


%% Delete file if exists
if exist(fn,'file')
    delete(fn)
end


%% Switch off XLS warning
ws = warning('off','MATLAB:xlswrite:AddSheet');


%% Get descriptives from function and write
descrtab = descriptives(data);
descrtab(:,[2 4]) = eachcell('str2num',descrtab(:,[2 4]));
xlswrite(fn,descrtab,'Descriptives','A1');


if ~isempty(output)
%% Get model fit stats and warnings from function and write
    [modfittab warns] = modelfittable(output,[],1);
    modfittab = [{'Model names','Deviance','Df','AICc','BIC','Runtime (sec)'};
        modfittab];
    xlswrite(fn,modfittab,'Model fits','A1');

%% Cycle through models
    for ctr=1:length(output)
        [nn np] = size(output(ctr).Minimum);
        if output(ctr).Options.NonparametricBootstrap>0
            minimum = output(ctr).NonparametricBootstrapMean;
            stderr = output(ctr).NonparametricBootstrapStdErr;
            partitle = 'Parameter estimates (nonparametric bootstrap)';
            stdtitle = 'Standard error estimates (nonparametric bootstrap)';
        elseif output(ctr).Options.ParametricBootstrap>0
            minimum = output(ctr).ParametricBootstrapMean;
            stderr = output(ctr).ParametricBootstrapStdErr;
            partitle = 'Parameter estimates (parametric bootstrap)';
            stdtitle = 'Standard error estimates (parametric bootstrap)';
        else
            minimum = output(ctr).Minimum;
            stderr = output(ctr).StdErr;
            partitle = 'Parameter estimates';
            stdtitle = 'Standard error estimates';
        end
        minimum = num2cell(minimum);
        stderr = num2cell(stderr);
        imags=~cellfun('isreal',stderr);
        stderr(imags)={'#IMAG'};
        parnames = {'Condition','a','Ter','eta','z','sZ','st','v'};
        if np>7
            parnames = [parnames {'pi' 'gamma'}];
        end
        shname = output(ctr).Name(1:min(end,31));

        % Write estimates
        xlswrite(fn,{output(ctr).Name},shname,'A1');
        xlswrite(fn,{partitle},shname,'A3');
        xlswrite(fn,parnames,shname,'A4');
        xlswrite(fn,(1:nn)',shname,'A5');
        xlswrite(fn,minimum,shname,'B5');
        xlswrite(fn,{stdtitle},shname,sprintf('A%u',6+nn));
        xlswrite(fn,parnames,shname,sprintf('A%u',7+nn));
        xlswrite(fn,(1:nn)',shname,sprintf('A%u',8+nn));
        xlswrite(fn,stderr,shname,sprintf('B%u',8+nn));

        % Write warnings
        oo = 0;
        if any(imags(:))
            xlswrite(fn,{['Some standard errors were not real numbers. This',...
                ' may indicate poor model fit.']},shname,sprintf('A%u',9+2*nn));
            oo=2;
        end

        xlswrite(fn,{'Warnings:'},shname,sprintf('A%u',9+oo+2*nn));
        xlswrite(fn,warns(ctr),shname,sprintf('B%u',9+oo+2*nn));
    end
end

%% Switch warning back to what it was
warning(ws.state,ws.identifier);


%% If called from base and on PC, launch Excel
if numel(dbstack)==1 && ispc
    try
        winopen(fn)
    catch
        warning('DMAT:outputtoxls:excelNotFound','Could not open MS Excel.')
    end
end