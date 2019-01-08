function varargout = modelfittable(output,fn,nums)
%MODELFITTABLE  Table with model properties
%   TABLE = MODELFITTABLE(OUTPUT), where OUTPUT is an array of DMAT output
%   structures, returns a cell array with some statistics about model fit.
%
%   MODELFITTABLE(OUTPUT) prints it to the screen.
%
%   MODELFITTABLE(OUTPUT,FILENAME), prints the table to the file FILENAME.
%
%   See also QPPLOT and EDFCDF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<1
    error('DMAT:modelfittable:notEnoughInputs',...
        'MODELFITTABLE requires at least one input variable.')
end

if nargin>2 && ~isempty(fn)
    if ispc
        fl = 'wt';
    else
        fl = 'w';
    end
    fid = fopen(fn,fl);
else
    fid = 1;
end

if nargin<3
    nums=false;
end

nmods = length(output);

models = cell(nmods,1);
devs = cell(nmods,1);
dfs = cell(nmods,1);
AICcs = cell(nmods,1);
LRs = cell(nmods,1);
BICs = cell(nmods,1);
warns = cell(nmods,1);

for ctr=1:nmods
    models(ctr)={output(ctr).Name};
    warns(ctr)={char(output(ctr).Warnings)};
    
    npb = output(ctr).Options.NonparametricBootstrap;
    ypb = output(ctr).Options.ParametricBootstrap;
    
    if npb>0
        fitvals = zeros(npb,1);
        aiccs = zeros(npb,1);
        lrats = zeros(npb,1);
        bics = zeros(npb,1);
        for iter = 1:npb
            fitvals(iter) = output(ctr).NonparametricBootstraps(iter).Fitvalue;
            aiccs(iter) = output(ctr).NonparametricBootstraps(iter).FitInfo.AICc;
            lrats(iter) = output(ctr).NonparametricBootstraps(iter).FitInfo.LogLikelihoodRatio;
            bics(iter) = output(ctr).NonparametricBootstraps(iter).FitInfo.BIC;
        end
        fitvalue = mean(fitvals);
        aicc = mean(aiccs);
        lr = mean(lrats);
        bic = sum(bics);
    elseif ypb>0        
        fitvals = zeros(ypb,1);
        aiccs = zeros(npb,1);
        lrats = zeros(npb,1);
        bics = zeros(ypb,1);
        for iter = 1:ypb
            fitvals(iter) = output(ctr).ParametricBootstraps(iter).Fitvalue;
            aiccs(iter) = output(ctr).NonparametricBootstraps(iter).FitInfo.AICc;
            lrats(iter) = output(ctr).NonparametricBootstraps(iter).FitInfo.LogLikelihoodRatio;
            bics(iter) = output(ctr).ParametricBootstraps(iter).FitInfo.BIC;
        end
        fitvalue = mean(fitvals);
        lr = mean(lrats);
        aicc = mean(aiccs);
        bic = sum(bics);
    else
        aicc = output(ctr).FitInfo.AICc;
        lr = output(ctr).FitInfo.LogLikelihoodRatio;
        fitvalue = output(ctr).Fitvalue;
        bic = output(ctr).FitInfo.BIC;
    end

    if ~nums
        devs(ctr)={num2str(fitvalue,'%9.4f')};
        dfs(ctr)={num2str(output(ctr).Df,'%3i')};
        AICcs(ctr)={num2str(aicc,'%9.4f')};
        LRs(ctr)={num2str(lr,'%9.4f')};
        BICs(ctr)={num2str(bic,'%9.4f')};
    else
        devs(ctr)={fitvalue};
        dfs(ctr)={output(ctr).Df};
        AICcs(ctr)={aicc};
        LRs(ctr)={lr};
        BICs(ctr)={bic};
    end
end


mfs=[models(:) devs(:) dfs(:) AICcs(:) BICs(:) LRs(:)];

if nargout>0
    varargout{1} = mfs;
    varargout{2} = warns(:);
else
    fprintf(fid,'\n     Model%17s',repmat(' ',1,17));
    fprintf(fid,'%13s','Deviance');
    fprintf(fid,'%5s','Df');
    fprintf(fid,'%13s','AICc');
    fprintf(fid,'%13s','BIC');
    fprintf(fid,'%13s\n','LR');
    for row = 1:size(mfs,1)
        fprintf(fid,'%3i  ',row);
        fprintf(fid,'%22s',parsetoline(mfs{row,1},22,1));
        fprintf(fid,'%13s',mfs{row,2});
        fprintf(fid,'%5s',mfs{row,3});
        fprintf(fid,'%13s',mfs{row,4});
        fprintf(fid,'%13s',mfs{row,5});
        fprintf(fid,'%13s\n',mfs{row,6});
    end
    fprintf(fid,'\n');
    if fid>2
        fclose(fid);
        disp('File written.')
    end
end