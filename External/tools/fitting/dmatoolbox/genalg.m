function [output,options]=genalg(options)
%GENALG  Generative algorithm for minimization
%   [OUTPUT,OPTIONS] = GENALG(OPTIONS), where OPTIONS is a full, valid
%   options structure for DMAT, returns a partial OUTPUT structure and an
%   amended OPTIONS structure.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.


%% Prepare options for first phase
nrops = optimset(optimset('fminunc'),...
    'FunValCheck','on',...
    'Display',options.Display,...
    'LargeScale','off',...
    'TolFun',10^(-options.ObjectiveDecimals-1),...
    'TolX',10^(-options.ParameterDecimals-1),...
    'MaxFunEvals',Inf,...
    'MaxIter',30000);
dform = sprintf('%%%i.%if',...
    options.ObjectiveDecimals+9,...
    options.ObjectiveDecimals+1);

%disp '      Starting generative algorithm'
objecfun = options.objecfun;

nmsops = optimset(optimset('fminsearch'),...
    'FunValCheck','on',...
    'Display',options.Display,...
    'LargeScale','off',...
    'TolFun',10^(-options.ObjectiveDecimals),...
    'TolX',10^(-options.ParameterDecimals),...
    'MaxFunEvals',Inf,...
    'MaxIter',250);

%% Begin generative minimization algorithm
output.Fitvalue=objecfun(options.controls.small);
fprintf(['\n      Guess     : ' dform '     (%s)\n'],output.Fitvalue,datestr(now))

%% Estimation stage with Nelder-Mead Simplex
nmsops.MaxIter = 200;
for rr=1:options.ShortSimplexRuns
    drawnow % for GUI
    [options.controls.small,output.Fitvalue,ign,output.Simplex(rr)] = ...
        fminsearch(objecfun,options.controls.small,nmsops); % Run simplex
    fprintf(['      Simplex %2u: ' dform '     (%s)\n'],rr,output.Fitvalue,datestr(now))
    %fprintf('          (Relative gain: %6.4f)\n',m-output.Fitvalue)
end
rr=options.ShortSimplexRuns; % in case this was 0

nmsops.MaxIter = options.MaxIter;
for rs=rr+1:(rr+options.LongSimplexRuns)
    drawnow % for GUI
    [options.controls.small,output.Fitvalue,ign,output.Simplex(rs)] = ...
        fminsearch(objecfun,options.controls.small,nmsops); % Run simplex
    fprintf(['      Simplex %2u: ' dform '     (%s)\n'],rs,output.Fitvalue,datestr(now))
    %fprintf('          (Relative gain: %6.4f)\n',m-output.Fitvalue)
end

%% Finalization
[options.controls.small,output.Fitvalue,options.errorflag(2),op,...
    gradient,output.Hessian] = fminunc(objecfun,options.controls.small,nrops);
fprintf(['      Final loss: ' dform '     (%s)\n\n'],output.Fitvalue,datestr(now))