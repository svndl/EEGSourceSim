function varargout = qtable(output,title)
% QTABLE  Shows concise summary of an output queue
%   QTABLE(OUTPUT), where OUTPUT is an array of output structures from
%   DMAT, prints a table with model deviances, degrees of freedom and
%   significance values of the sequential changes in model fit, as well as
%   AICc (Small Sample AIC) and BIC values.
%
%   QTABLE(OUTPUT,TITLE) adds a title to the table.
%
%   TABLE = QTABLE(OUTPUT) prints nothing but returns the same table in a
%   matrix format.
%
%   See also RUNQUEUE and DMATGUI.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Fixed typo in AICc computation.
%     Was:    AICc = -2*Lm + 2*k + 2.*k.*(k-1)./(n-k-1);
%     Is now: AICc = -2*Lm + 2*k + 2.*k.*(k+1)./(n-k-1); (from ABSFIT)

if ~isstruct(output)
    if ischar(output)
        if exist(output,'file')
            c = load(output);
            fl = fieldnames(c);
            output = c.(fl{1});
        else
            error('DMAT:qtable:badOutputStructure',['%s is not a file '...
                'on the MATLAB path.'],output)
        end
    else
        error('DMAT:qtable:badOutputStructure',['Argument OUTPUT '...
            'should be an output structure.'])
    end
end

nm = numel(output);

for a = 1:nm
    dev(a) = output(a).Fitvalue;
    df(a) = output(a).Df;
    AI(a) = output(a).FitInfo.AICc;
    BI(a) = output(a).FitInfo.BIC;
end
ddev = [NaN -diff(dev)];
ddf = [NaN diff(df)];
sig = chi2test(ddev,ddf);


if nargout
    varargout{1} = [dev' ddev' df' ddf' AI' BI' sig'];
else
    line = ['   ' repmat('-',1,71) '\n'];
    tl = line;
    if nargin==2
        tl(7:length(title)+8) = [' ' title ' '];
    end
    fprintf(tl)
    fprintf('%14s %11s %4s %6s %11s %11s %9s\n','Deviance',...
        'd-Deviance','df','d-df','AICc','BIC','p')
    fprintf(line)
    fprintf('%14.4f %11.4f %4i %6i %11.4f %11.4f %9.5f\n',...
        [dev' ddev' df' ddf' AI' BI' sig']')
    fprintf(line)
end

%{
%% Obsolete code
function varargout = qtable_n(output,title)
% QTABLE  Shows concise summary of an output queue
%   QTABLE(OUTPUT), where OUTPUT is an array of output structures from
%   DMAT, prints a table with model deviances, degrees of freedom and
%   significance values of the sequential changes in model fit, as well as
%   AICc (Small Sample AIC) and BIC values.
%
%   QTABLE(OUTPUT,TITLE) adds a title to the table.
%
%   TABLE = QTABLE(OUTPUT) prints nothing but returns the same table in a
%   matrix format.
%
%   See also RUNQUEUE and DMATGUI.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if ~isstruct(output)
    if ischar(output)
        if exist(output,'file')
            c = load(output);
            fl = fieldnames(c);
            output = c.(fl{1});
        else
            error('DMAT:qtable:badOutputStructure',['%s is not a file '...
                'on the MATLAB path.'],output)
        end
    else
        error('DMAT:qtable:badOutputStructure',['Argument OUTPUT '...
            'should be an output structure.'])
    end
end

dev = [output.Fitvalue];
ddev = [NaN -diff(dev)];
df = [output.Df];
ddf = [NaN diff(df)];
sig = chi2test(ddev,ddf);
tmp = [output.OutlierReport];
n = sum([tmp.use]);
AI = aic(dev,df,n);
BI = bic(dev,df,n);
if nargout
    varargout{1} = [dev' ddev' df' ddf' AI' BI' sig'];
else
    line = ['   ' repmat('-',1,71) '\n'];
    tl = line;
    if nargin==2
        tl(7:length(title)+8) = [' ' title ' '];
    end
    fprintf(tl)
    fprintf('%14s %11s %4s %6s %11s %11s %9s\n','Deviance',...
        'd-Deviance','df','d-df','AICc','BIC','p')
    fprintf(line)
    fprintf('%14.4f %11.4f %4i %6i %11.4f %11.4f %9.5f\n',...
        [dev' ddev' df' ddf' AI' BI' sig']')
    fprintf(line)
end

function a = aic(dev,k,n), a = dev + 2*k + 2.*k.*(k-1)./(n-k-1);
function a = bic(dev,k,n), a = dev + k.*log(n);
%}