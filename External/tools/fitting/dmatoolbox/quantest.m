function output=quantest(varargin)
%QUANTEST  Fits the diffusion model to a data represented by quantiles only
%   OUTPUT = QUANTEST(DSTR,OPTIONS), where DSTR is a structure with all the
%   fields:
%            .CorrectObs          .Number of observations in each X=1 bin
%            .IncorrectObs        .Number of observations in each X=0 bin
%            .CorrectEdges        .Edges of X=1 time bins (in sec)
%            .IncorrectEdges      .Edges of X=0 time bins (in sec)
%            .Lower               .Minimal RT (optional)
%            .Upper               .Maximal RT (optional)
%   works just like MULTIESTV4, with the exception of outlier analyses.
%   Mixture model outlier analysis requires that you provide the .Lower and
%   .Upper fields. EWMA is not possible.
%   QUANTEST cannot make automatic use of the EZDIFF algorithm, if you can
%   provide EZDIFF results as starting points for the parameter estimation,
%   you should do so explicitly by providing the .Guess field in the
%   OPTIONS structure.
%
%   See also MULTIESTV4, RUNQUEUE_Q.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Forward compatibility notice: This function will likely be absorbed into
%  MULTIESTV4 in a future DMAT release.

%% Define default options structure
defopts = struct(...
    'DesignMatrix',{repmat({'1'},1,7)},...
    'Display','off',...
    'EWMA',struct('l',.01,'L',1.5,'s',.5),...
    'EstimationMethodScalar',5,...    'Fastdm',0,...
    'FixedBinEdges',[0.30 0.36 0.42 0.52 0.80;0.38 0.47 0.56 0.70 1.00],...
    'FixedValues',[],...
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
elseif nargin==1 % If none given, use default options structure
    data = varargin{1};
    options = defopts;
elseif nargin==2
    data = varargin{1};
    options = varargin{2};
elseif nargin>2
    data = varargin{1};
    options = pairs2struct(varargin{2:end});
end

%% Check if user submitted a queue
if length(options)>1
    output = runqueue_q(data,options);
    return
end

%% Start timer
startingtime = clock;

%% Check if input is valid
options.nsets = numel(data);
options = inpcheck(options);

%% Prepare for minimization
% Get formatted data from input structure
[dcel dcel2] = str2dcel(data);

% Generate guess
options = generateguess(options,dcel);
options.controls.large = options.Guess;

% Censor unnecessary parameters
options.controls = smaller(options.controls);
options.controls.display = options.Display;

% Define the objective function
options.objecfun = @(x) (multifitv4(x,dcel2,options.controls,options.method));
%options.objecfun = @(x) (multifitv7(x,splitdata(data,1),options.controls));

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
        options.Guess(:,5) = .95*min(output.Minimum(:,4),...
            output.Minimum(:,1)-output.Minimum(:,4));
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
        options.Guess(:,6) = options.Guess(:,2)/2;
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
% parameter set (mainly for simulation analyses)
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
output.Options = options;
output.Df = length(options.controls.small);
output.DesignVector = options.controls.small;
output.FitInfo = absfit(output.Minimum,dcel2,output.Df);
output.Options.controls.dcel2 = dcel2;
output.OutlierReport.use = true(sum(cellfun(@(x)x{7},dcel2)),1);

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

%%
function [badter badsz badeta badst]= crediblesolution(minimum,keepers)
% Is the solution at the edge of the parameter space?
badter = (minimum(:,2)<.02 & keepers(:,2));
badsz = (minimum(:,5)<minimum(:,4)/100) & keepers(:,5) & keepers(:,4);
badeta = (minimum(:,3)<.0005) & keepers(:,3) & keepers(:,7);
badst = (minimum(:,6)<minimum(:,2)/100) & keepers(:,6) & keepers(:,2);