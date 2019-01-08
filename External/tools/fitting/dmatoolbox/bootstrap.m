function [op,minimum] = bootstrap(data,output,iter,parametric)
%BOOTSTRAP  Parametric or nonparametric bootstrap of DMAT results
%   [OP,MINIMUM] = BOOTSTRAP(DATA,OUTPUT,ITER,PARAMETRIC), where DATA is a
%   valid data set and OUTPUT is either (1) a DMAT output structure or (2)
%   a structure with a DMAT options structure in one field called .Options
%   and an initial guess in a field called .Minimum. ITER is the number of
%   iterations and PARAMETRIC is a logical value indicating if you want
%   parametric bootstrap (true for parametric, false for nonparametric).
%
%   All input arguments except DATA are optional. Default settings are 20
%   iterations of nonparametric bootstrap runs, with the default options
%   from MULTIESTV4.
%   If you specify parametric bootstraps, the OUTPUT variable is required
%   and cannot be empty.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<4 || isempty(parametric)
    parametric=0;
elseif ~isscalar(parametric)
    parametric = parametric~=0;
end

if nargin<3 || isempty(iter)
    iter=20;
elseif ~isscalar(iter) || iter<0
    error('DMAT:bootstrap:badInput',...
        'Third input variable to BOOTSTRAP should a positive scalar.');
end

if nargin<2 || isempty(output)
    if parametric
        error('DMAT:bootstrap:parametricBootstrapRequiresGuess',...
            'For parametric bootstraps, you need to input parameters.');
    else
        output = struct('Options',multiestv4,'Minimum',[]);
    end
elseif ~isstruct(output) || ~isfield(output,'Options') || ~isfield(output,'Minimum')
    error('DMAT:bootstrap:badInput',...
        'Second input variable to BOOTSTRAP should a DMAT output structure.');
end

if ~nargin
    error('DMAT:bootstrap:notEnoughInputs',...
        'BOOTSTRAP requires at least one input variable.')
end

options = output.Options;

options.Guess = output.Minimum;
options.InBootstrap = true;

ncond = options.nsets;
dcel = splitdata(data);
ns = cellfun('length',dcel);

if parametric
    for it=1:iter
        e = multisimul(output.Minimum,ns);
        fprintf('\n     Parametric bootstrap iteration %i\n',it);
        op(it) = multiestv4(e,options); %#ok
        minimum(:,:,it) = op(it).Minimum;
    end
else
    subsample = zeros(sum(ns),3);
    ns = [0;ns(:)]';
    for it=1:iter
        for cnd = 1:ncond
            ind = ceil(rand(ns(cnd+1),1)*ns(cnd+1));
            subsample(sum(ns(1:cnd))+1:sum(ns(1:(cnd+1))),2:3) = dcel{cnd}(ind,:);
            subsample(sum(ns(1:cnd))+1:sum(ns(1:(cnd+1))),1) = cnd;
        end
        fprintf('\n     Nonparametric bootstrap iteration %u\n',it);
        temp = multiestv4(subsample,options);
        op(it) = temp; %#ok
        minimum(:,:,it) = op(it).Minimum;
    end
end