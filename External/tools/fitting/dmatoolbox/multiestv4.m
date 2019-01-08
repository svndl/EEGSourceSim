function output=multiestv4(varargin)
%MULTIESTV4  Fits the diffusion model to a data set
%   OUTPUT = MULTIESTV4(DATA,OPTIONS), where DATA is a N-by-3 data set with
%   columns [condition response seconds] and OPTIONS is an (optional)
%   structure with model and optimization information.
%   Note: If OPTIONS is an array of structures, the function RUNQUEUE will
%   be invoked.
%
%   OUTPUT = MULTIESTV4(DATA) uses a default options structure.
%
%   MULTIESTV4 also works with property-value pairs. For example:
%   OUTPUT = MULTIESTV4(DATA,'DesignMatrix',repmat({orthpoly(4,0)},1,7))
%   will use the specified design matrices and use default options for the
%   rest.
%
%   Special uses:
%     DEFOPTS = MULTIESTV4(N) returns a vector with N default options
%     structures.
%     OPTIONS = MULTIESTV4(DATA,OPTIONS), where the .NoFitting field in
%     OPTIONS is set to true, doesn't estimate parameters, but returns an
%     OPTIONS structure with the correct objective function in the
%     .objecfun field.
%
%   For further information regarding the use of this function, please see
%   the DMAT manual. For technical information or an introduction to
%   diffusion model analysis, please see the website, or type DMATREF for
%   the references to the software and accompanying article.
%
%   If you did not obtain this m-file by downloading it from the DMAT
%   website, please find the most recent version there and download it and
%   relevant information for free.
%
%   See also DMATSITE and RUNQUEUE.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Scalar input now returns array of options structures.
%  Edit 0.4: Stored options structure now contains dcel and dcel2.
%  Edit 0.5: Added option to fit B, not Z.
%  Edit 0.5: Changed jump strategy because it would sometimes jump out of
%            the parameter space.

%% Define default options structure
defopts = struct(...
    'DesignMatrix',{repmat({'1'},1,7)},...
    'Display','off',...
    'EWMA',struct('l',.01,'L',1.5,'s',.5),...
    'EstimationMethodScalar',5,...    'Fastdm',0,...
    'FixedBinEdges',[0.30 0.36 0.42 0.52 0.80;0.38 0.47 0.56 0.70 1.00],...
    'FixedValues',[],...
    'FitBnotZ',0,...
    'Guess',[],...
    'GuessMethodScalar',1,...
    'LongSimplexRuns',1,...
    'MaxIter',5000,...    'ModelInfo',0,...
    'Name','No name given',...
    'NoFitting',0,...
    'NonparametricBootstrap',0,...
    'ObjectiveDecimals',7,...
    'OutlierMax',[],...
    'OutlierMin',[],...
    'OutlierTreatment','None',...
    'ParameterDecimals',7,...
    'ParametricBootstrap',0,...
    'Percentiles',[10 30 50 70 90;10 30 50 70 90],...
    'ShortSimplexRuns',3,...
    'SpecificBias',[]);

if ~nargin % If no input arguments given, return default options structure
    output = defopts;
    return
elseif nargin==1
    % Or return a vector of them
    if isscalar(varargin{1})
        output = repmat(defopts,varargin{1},1);
        return
    end
    % If no options structure given, use default options structure
    data = varargin{1};
    options = defopts;
elseif nargin==2
    data = varargin{1};
    options = varargin{2};
elseif nargin>2
    data = varargin{1};
    options = pairs2struct(varargin{2:end});
end

% if numel(options)==1 && isfield(options,'Fastdm') && options.Fastdm
%     type('fdmgnu.txt')
%     output = fastdm(data);
%     return
% end

%% Check if user submitted a queue
if length(options)>1
    output = runqueue(data,options);
    return
end

%% Start timer
startingtime = clock;

%% First, get number of conditions from data
switch size(data,2)
    case 2
        if isvaliddataset(data,2)
            options.nsets = 1;
            data = [ones(size(data,1),1) data];
        else
            error('DMAT:multiestv4:invalidDataSet',...
                'Data should be in format [response,RT] or [condition,response,RT]')
        end
    case 3
        if isvaliddataset(data,3)
            options.nsets = length(unique(data(:,1)));
        else
            error('DMAT:multiestv4:invalidDataSet',...
                'Data should be in format [response,RT] or [condition,response,RT]')
        end
    otherwise
        error('DMAT:multiestv4:invalidDataSet',...
            'Data should be in format [response,RT] or [condition,response,RT]')
end

%% Check if input is valid
options = inpcheck(options);

%% Treat outliers, generate report
[data,oltrreport] = outliertreatment(data,options);

%% Split data into cells
% This function splits one dataset into a dataset for each condition.
dcel = splitdata(data);

%% Prepare for minimization
% This function takes a split dataset and calculates the necessary
% (quantiles etc) statistics on each condition
dcel2 = processdata(dcel,options.Percentiles,options.FixedBinEdges);

% Generate guess
options = generateguess(options,dcel);
options.controls.large = options.Guess;
options.controls.fitbnotz = options.FitBnotZ;

% Censor unnecessary parameters
options.controls = smaller(options.controls);
options.controls.display = options.Display;

% Define the objective function
options.objecfun = @(x) (multifitv4(x,dcel2,options.controls,options.method));
%options.objecfun = @(x) (multifitv7(x,splitdata(data,1),options.controls));

%% If NoFitting was on, stop here
if options.NoFitting
    options.OutlierReport = oltrreport;
    output = options;
    return
end

%% Run generative minimization algorithm
[output options]=genalg(options);

%% Rebuild full parameter matrix
options.controls = bigger(options.controls);
output.Minimum = options.controls.large;

%% Check solution
[badter badsz badeta badst]= crediblesolution(output.Minimum,...
    options.controls.keepers);
options.retrys=0;

bestoutput = output;

while any(any([badter badsz badeta badst])) && options.retrys<1
    options.retrys=options.retrys+1;
    options.Guess = output.Minimum;
    
    if any(badsz)
        fprintf('      The recovered sZ parameters are suspect.\n')
        options.Guess(:,5) = .95*min(min([output.Minimum(:,4),...
            output.Minimum(:,1)-output.Minimum(:,4)]));
    end
    if any(badeta)
        fprintf('      The recovered eta parameters are suspect.\n')
        options.Guess(:,3) = .2;
    end
    if any(badter)
        fprintf('      The recovered Ter parameters are suspect.\n')
        options.Guess(:,2) = .25;
    end
    if any(badst)
        fprintf('      The recovered st parameters are suspect.\n')
        options.Guess(:,6) = .9*min(options.Guess(:,2));
    end
    fprintf('      Trying again.\n')
    options.controls.large = options.Guess;
    options.controls = smaller(options.controls);

    %% Restart generative minimization algorithm
    [output options]=genalg(options);
    options.controls = bigger(options.controls);
    output.Minimum = options.controls.large;
    
    [badter badsz badeta badst]= crediblesolution(output.Minimum,...
        options.controls.keepers);
    
    if output.Fitvalue<bestoutput.Fitvalue
        bestoutput = output;
    end
end

if any(any([badter badsz badeta badst]))
    fprintf('      The last convergence point was still a suspect result.\n')
end
if any(bestoutput.Minimum~=output.Minimum)
    fprintf('      Returning to the best point found and giving up.\n\n')
end
output = bestoutput;

%% Process Hessian
[output.StdErr warn]= processhessian(output.Hessian,options);


%% Finish output
% If an extra parameter set was passed, calculate badness-of-fit for that
% parameter set (mainly for simulation analyses). Doesn't work particularly
% well.
if isfield(options,'TruePar') && ~isempty(options.TruePar)
    output.Badness = 0;
    for ctr=1:options.nsets
        output.Badness = output.Badness + ...
            fitdiffv13(options.TruePar(ctr,:),dcel2{ctr},options.method);
    end
end
if ~options.errorflag
    warn{end+1} = 'Simplex optimization did not converge.';
end
if ~options.errorflag(2)
    warn{end+1} = 'Newton-Raphson optimization did not converge, solution may not be a minimum.';
end
if isempty(warn)
    warn = {'No warnings.'};
end
output.Warnings = warn;
output.Name = options.Name;
output.Time = etime(clock,startingtime);
output.OutlierReport = oltrreport;

for ctr=1:options.nsets
    if ~isempty(dcel2{ctr}{1})
        options.yQ(ctr,1:length(dcel2{ctr}{1})) = dcel2{ctr}{1};
    else
        options.yQ(ctr,1:length(dcel2{ctr}{1})) = 0;
    end
    if ~isempty(dcel2{ctr}{2})
        options.nQ(ctr,1:length(dcel2{ctr}{2})) = dcel2{ctr}{2};
    else
        options.nQ(ctr,1:length(dcel2{ctr}{2})) = 0;
    end
end
options.dcel = dcel;
options.dcel2 = dcel2;
output.Options = options;
output.Df = length(options.controls.small);
output.DesignVector = options.controls.small;
output.FitInfo = absfit(output.Minimum,dcel2,output.Df);
output.Options.controls.dcel2 = dcel2;

%% Optionally, bootstrap confidence intervals
output.ParametricBootstrapStdErr =[];
output.NonparametricBootstrapStdErr =[];
output.ParametricBootstrapMean=[];
output.ParametricBootstraps=[];
output.NonparametricBootstrapMean=[];
output.NonparametricBootstraps=[];

if ~options.InBootstrap
    if options.ParametricBootstrap>0 % Parametric
        [output.ParametricBootstraps minim] = ...
            bootstrap(data,output,options.ParametricBootstrap,true);
        output.ParametricBootstrapStdErr = std(minim,0,3);
        output.ParametricBootstrapMean = mean(minim,3);
    elseif options.NonparametricBootstrap>0 % Nonparametric
        [output.NonparametricBootstraps minim] = ...
            bootstrap(data,output,options.NonparametricBootstrap,false);
        output.NonparametricBootstrapStdErr = std(minim,0,3);
        output.NonparametricBootstrapMean = mean(minim,3);
    end
end

output = orderfields(output);


function [badter badsz badeta badst]= crediblesolution(minimum,keepers)
% Is the solution at the edge of the parameter space?
badter = (minimum(:,2)<.02 & keepers(:,2));
badsz = (minimum(:,5)<minimum(:,4)/100) & keepers(:,5) & keepers(:,4);
badeta = (minimum(:,3)<.0005) & keepers(:,3) & keepers(:,7);
badst = (minimum(:,6)<minimum(:,2)/100) & keepers(:,6) & keepers(:,2);