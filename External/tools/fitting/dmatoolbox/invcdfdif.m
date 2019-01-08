function rts = invcdfdif(perc,rep,output,precision)
%INVCDFDIF  Inverse of the diffusion defective cumulative density function
%  X = INVCDFDIF(PERC,REP,OUTPUT,PRECISION), for output structure OUTPUT
%  and response REP, returns the response times associated with the
%  percentiles in PERC (for each condition), to a precision of PRECISION.
%  The values in PERC should be for the defective CDF (that is, between 0
%  and 1).
%
%  X = INVCDFDIF(PERC,REP,OUTPUT), uses a precision of 1e-4.
%
%  X = INVCDFDIF(PERC,REP,PARAM), for parameter vector PARAM and response
%  REP, returns the response times associated with the percentiles in PERC.
%
%  See also CDFDIF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command.
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.4: Rewrote from scratch. Now uses recursive lookup table method
%   instead of (unreliable and inefficient) optimization method. Old code
%   is documented below.

if nargin<4
    precision = 1e-4;
end
if nargin<3
    error('DMAT:invcdfdif:notEnoughInputs',...
        'INVCDFDIF requires at least three input arguments.');
end
if ~isscalar(rep) || ~(rep==0 || rep==1)
    error('DMAT:invcdfdif:BadResponse','Response REP must be 0 or 1.')
end
if ~numel(perc)
    error('DMAT:invcdfdif:NoPercentile',...
        'Percentile vector PERC should not be empty.')
end
if any(perc<0)
    warning('DMAT:invcdfdif:negativeProbs',...
        'Negative probabilities detected in INVCDFDIF input.')
end

if isnumeric(output)
    havestruct = false;
    param = output;
else
    if numel(output)>1
        error('DMAT:invcdfdif:onlyOneStructureAllowed',...
            'You can only input a single output structure to INVCDFDIF.')
    end
    havestruct = true;
    param = output.Minimum;
end
if ~isgood(param) || (~havestruct && ~size(param,1)==1)
%    error('DMAT:invcdfdif:badParameterVector',...
%        'Input variable PARAM does not seem to be a good parameter vector.')
    warning('DMAT:invcdfdif:badParameterVector',...
        'Input variable PARAM does not seem to be a good parameter vector.')
    rts = nan(size(perc));
    return
end

% Prepare output matrix
ncond = size(param,1);
nperc = size(perc,2);
rts = zeros(ncond,nperc); 
perc = repmat(perc,ncond,1);

% Loop over inputs
for ctr = 1:ncond
    [c px]=cdfdif(0,rep,param(ctr,:)); %#ok
    px = 2*rep*px-rep-px+1; % Get right edge
    rts(ctr,perc(ctr,:)<=0) = NaN; % Left of zero = NaN
    rts(ctr,1<=perc(ctr,:)) = Inf; % Right of domain = Inf
    tofind = rts(ctr,:)==0; % Identify the rest
    rts(ctr,tofind) = subivc(perc(ctr,tofind),param(ctr,:),rep,...
        param(2)-param(6)/2-.001,5,px,precision); % Look-up table for perc
end

%% Subfunction SUBIVC
% Finds RT within a range that best approaches cdfidf^{-1}(perc).
% PER are the percentiles, PAR the parameter vector, REP the response, MI
% the lower bound, MA the upper bound, PX the conditional probability of
% REP and PRECISION the resolution of the inverse.
function r = subivc(per,par,rep,mi,ma,px,precision)
nsp = 12; % Number of points in the domain to sample
domain = linspace(mi,ma,nsp); % Points in the domain
table = cdfdif(domain,rep,par)./px; % Their associated CDF values
nperc = length(per); % Number of percentiles to approach
r = zeros(1,nperc); % Output vector

for cc=1:nperc 
    dd = abs(per(cc)-table); % Difference from lookup points
    [match ind] = min(dd); %#ok % Location of smallest difference
    lob = domain(max(1,ind-1)) ; % One location lower
    upb = domain(min(ind+1,nsp)); % One location higher
    if abs(upb-lob)>2*precision % If resulting resolution too crude
        r(cc) = subivc(per(cc),par,rep,... % Sample from smaller domain
            lob,upb,px,precision); % With same settings
    else
        r(cc) = domain(ind); % Else return recovered value
    end
end


%% This is the old (now retired) code
%{
function rts = invcdfdif(perc,rep,param)
%INVCDFDIF  Inverse of the diffusion defective cumulative density function
%  X = INVCDFDIF(PERC,REP,PARAM), for parameter vector PARAM and response
%  REP, returns the response times associated with the percentiles in PERC.
%  The values in PERC should be possible for the combination of parameters
%  and response (i.e., should be larger than 0 and smaller than the
%  marginal probability of the response REP).
%
%  See also CDFDIF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

% Obsolete function
warning('DMAT:invcdfdif:obsolete',...
    'Function INVCDFDIF is obsolete. Use CDFDIFINV instead.')
    
% Check input
if nargin<3
    error('DMAT:invcdfdif:notEnoughInputs',...
          'INVCDFDIF requires three input arguments.');
end
if ~isscalar(rep) || ~(rep==0 || rep==1)
    error('DMAT:invcdfdif:BadResponse','Response REP must be 0 or 1.')
end
if ~numel(perc)
    error('DMAT:invcdfdif:NoPercentile',...
        'Percentile vector PERC should not be empty.')
end
if ~size(param,1)==1 || ~isgood(param)
    error('DMAT:invcdfdif:badParameterVector',...
        'Input variable PARAM does not seem to be a good parameter vector.')
end

if any(perc<0)
    warning('DMAT:invcdfdif:negativeProbs',...
        'Negative probabilities detected in INVCDFDIF input.')
end

% Treat edge and out-of-bounds cases
rts = repmat(-9999,size(perc)); % Prepare output vector with known flags
[c px]=cdfdif(0,rep,param); px = 2*rep*px-rep-px+1; %#ok % Get right edge
rts(perc<0) = NaN; % Left of zero = NaN
rts(px<=perc) = Inf; % Right of domain = Inf
tofind = rts==-9999; % Identify the rest, where flags are untouched
perc = perc(tofind);

[r,v,exitflag] = fminunc(@objective,log(param(2)+perc),...
    struct('LargeScale','off',...
    'Display','off',...
    'TolX',1e-12,...
    'TolFun',1e-16,...
    'InitialHessType','identity')); %#ok
if ~exitflag
    beep
    disp 'Switching to NMS!'
    pause(1)
    r = fminsearch(@objective,r,...
        struct('Display','off',...
        'TolX',1e-12,...
        'TolFun',1e-16));
end

r = param(2)-param(6)/2+.001+exp(r);
rts(tofind)=r;

    function l=objective(logrts)
        r_rts = param(2)-param(6)/2+.001+exp(logrts);
        r_percs = cdfdif(r_rts,rep,param);
        l=r_percs-perc;
        l=l*l';
    end

end
%}