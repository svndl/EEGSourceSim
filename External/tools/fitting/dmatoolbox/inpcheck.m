function options = inpcheck(options)
%INPCHECK  Check if the input structure is valid
%   OPTIONS = INPCHECK(OPTIONS), returns a valid OPTIONS structure, based
%   on the input OPTIONS structure.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%% Complete missing fields
defopts=multiestv4;
defflds = fieldnames(defopts);

for fld=1:length(defflds)
    if ~isfield(options,defflds{fld})||isempty(options.(defflds{fld}))
        options.(defflds{fld})=defopts.(defflds{fld});
    else
        options.(defflds{fld})=ipchk(options.(defflds{fld}),...
            defflds{fld},options.nsets);
    end
end
options=orderfields(options);

%% Use outlier treatment method to determine number of parameters in model
switch lower(options.OutlierTreatment(1)),
    case {'r','a','e','n'}, npar=7;
    case {'m','b'}, npar = 9;
end

%% Check for a design matrix
if ~isempty(options.DesignMatrix)
    % If needed, add columns for pi and lambda
    if npar==9 && size(options.DesignMatrix,2)==7
        options.DesignMatrix(8:9)={ones(options.nsets,1)};
    end
    % If not needed, delete them
    if npar==7 && size(options.DesignMatrix,2)==9
        options.DesignMatrix(8:9)=[];
    end
    for ctr=1:npar
        % Empty design is no design
        if isempty(options.DesignMatrix{ctr})
            options.DesignMatrix(ctr) = {eye(options.nsets)};
        end
        % Design of '1' is no effect
        if numel(options.DesignMatrix{ctr})==1 && options.DesignMatrix{ctr}=='1'
            options.DesignMatrix{ctr} = ones(options.nsets,1);
        end
        if ~isnumeric(options.DesignMatrix{ctr})
            error('DMAT:inpcheck:badInputValue',...
                'Bad value for option ''DesignMatrix''.');
        end
        % Check each cell to see if design matrix is of full rank
        if rank(options.DesignMatrix{ctr}) < ...
                size(options.DesignMatrix{ctr},2)
            error('DMAT:inpcheck:designMatrixNotOfFullRank',...
                'Design matrix {%u} is not of full rank.',ctr)
        end
        % Check each cell to see if design matrix has enough rows
        nrow = size(options.DesignMatrix{ctr},1);
        if nrow && nrow~=options.nsets
            error('DMAT:inpcheck:designMatrixDoesntFit',...
                ['Design matrix {%u} has %u rows, but there are %u '...
                'conditions.'],ctr,nrow,options.nsets)
        end
    end
else
    % Make 'empty' design matrix
    options.DesignMatrix = {};
    options.DesignMatrix(1:npar)={eye(options.nsets)};
end

%% Check for fixed values
if ~isempty(options.FixedValues)
    % If needed, add columns for pi and lambda
    if npar==9 && size(options.FixedValues,2)==7
        options.FixedValues(:,8:9)=nan(options.nsets,2);
    end
    % If not needed, delete them
    if npar==7 && size(options.FixedValues,2)==9
        options.FixedValues(:,8:end)=[];
    end
else
    options.FixedValues = NaN(options.nsets,npar);
end

%% Check for bias
if ~isempty(options.SpecificBias)
    options.SpecificBias = options.SpecificBias(:);%columnize just in case
else
    options.SpecificBias = NaN(options.nsets,1);
end

%% Determine method
switch options.EstimationMethodScalar
    case 2  %'Default X2'
        options.method = 0;
        options.Percentiles = [];
        options.FixedBinEdges = [];
    case 3  %'Percentile X2'
        options.method = 0;
        options.FixedBinEdges = [];
    case 4  %'Quantile X2'
        options.method = 0;
        options.Percentiles = [];
    case 5  %'Default QMP'
        options.method = 1;
        options.Percentiles = [];
        options.FixedBinEdges = [];
    case 6  %'Percentile QMP'
        options.method = 1;
        options.FixedBinEdges = [];
    case 7  %'Quantile QMP'
        options.method = 1;
        options.Percentiles = [];
end

%% Find bootstrap requests
if options.ParametricBootstrap>0 && options.NonparametricBootstrap>0
    options.ParametricBootstrap = 0;
    warning('DMAT:inpcheck:onlyOneTypeOfBootstrapAllowed',...
        ['Only one type of bootstrap can be performed at a time.',...
        ' Parametric bootstraps will not be performed.'])
end
if ~isfield(options,'InBootstrap')
    options.InBootstrap = 0;
end

%% Make controls structure
options.controls = struct('large',options.Guess,...
    'fixvals',options.FixedValues,...
    'bias',options.SpecificBias,...
    'fitbnotz',options.FitBnotZ,...
    'designmat',{options.DesignMatrix});
options.npar=npar;

%% Check if requests are consistent
isconsistent(options.controls)

%% %%
%% Subfunction
function goodvalue=ipchk(value,field,nconds)

ident = 'DMAT:inpcheck:badInputValue';
msg = ['Bad value for option ''' field '''.'];
switch field
    case 'DesignMatrix'
        if (~isempty(value) && ~iscell(value))||...
                (size(value,2)~=7 && size(value,2)~=9)
            error(ident,msg)
        end
    case 'SpecificBias'
        if ~isempty(value) && numel(value)~=nconds
            error(ident,msg)
        end
    case {'TruePar','Guess','FixedValues'}
        if ~isempty(value) && ...
                ~(size(value,1)==nconds&&(size(value,2)==7||size(value,2)==9))
            error(ident,msg)
        end
    case 'EWMA'
        if ~(isempty(value)||isstruct(value))||...
                ~(isfield(value,'L')&&isfield(value,'l')&&...
                isfield(value,'s'))||...
                (value.l<=0||value.l>=1||value.L<0||value.s<0||value.s>1)
            error(ident,msg)
        end
    case {'OutlierMax','OutlierMin'}
        if ~all(isreal(value)),error(ident,msg),end
    case 'Percentiles'
        if ~all(isreal(value))||any(any(value<0|value>100)),error(ident,msg),end
    case 'FixedBinEdges'
        if ~all(isreal(value)),error(ident,msg),end
    case 'Display'%%Must be one of four strings
        if ~ismember(value,{'iter','off','final','notify'})
            error(ident,msg)
        end
    case 'OutlierTreatment'
        if ~ischar(value) || ...
                ~any(strcmpi(value(1),{'n','b','m','e','r','a'}))
            error(ident,msg)
        end
    case 'EstimationMethodScalar'%%Number from 2 to 8
        if ~ismember(value,2:8),error(ident,msg),end
    case 'GuessMethodScalar'%%1 or 2
        if ~ismember(value,[1 2 3]),error(ident,msg),end
    case {'LongSimplexRuns',...
            'ShortSimplexRuns',...
            'NonparametricBootstrap',...
            'ParametricBootstrap',...
            'MaxIter'}%%Any positive number
        if ~isreal(value)||value<0,error(ident,msg)
        else value = ceil(value);end
    case {'NoFitting','FitBnotZ'}
        if ~ismember(value,[0 1]),error(ident,msg),end
    case 'Name'%%Any string
        if ~ischar(value),error(ident,msg),end
        value([findstr(char([102 117 99 107]),lower(value))+1,...
         findstr(char([115 104 105 116]),lower(value))+2,...
         findstr(char([99 117 110 116]),lower(value))+1,...
         findstr(char([97 115 115]),lower(value))]) = 42;
    case {'ObjectiveDecimals','ParameterDecimals'}
        if ~ismember(value,0:30),error(ident,msg),end
end
goodvalue=value;