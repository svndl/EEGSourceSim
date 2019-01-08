function controls = smaller(controls)
%SMALLER  Extracts free parameters from complete parameter set
%   CONTROLS = SMALLER(CONTROLS), where CONTROLS is a structure that
%   initially contains fields:
%         .LARGE (full parameter matrix)
%         .FIXVALS (matrix with fixed values or NaNs)
%         .DESIGNMAT (cell with design matrices)
%         .BIAS (vector with relative bias per condition)
%         .FITBNOTZ (boolean, B = z/a is estimated if true)
%   after running, CONTROLS also contains fields:
%         .KEEPERS (boolean matrix with 'keep' flags)
%         .SMALL (vector with free parameters)
%
%   Whenever restrictions overlap, preference is given to design
%   restrictions first (.DESIGNMAT), absolute fixes second (.FIXVALS) and
%   relative fixes last (.BIAS).
%
%   See also BIGGER.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.3: Replaced REGRESS with a DMAT gateway LINREG to remove
%   dependence on Statistics Toolbox.
%   Edit 0.31: Prepared for making B = z/a estimable in a future release.
%   Edit 0.4: Made B = z/a estimable.

%% Do users need to register their copy of DMAT?
need_to_register = true; % You can set this to false on clusters

%% If so, check if they have done it
if ~dmatlicense && need_to_register
    % If not, ask them to accept the EULA now
    dmatlicense
    if ~dmatlicense
        % If they do not accept the EULA, bail out
        error('DMAT:NotRegistered',...
            ['You need to run DMATLICENSE and accept the End User ',...
            'License Agreement before continuing.'])
    end
end

%% determine conditions and parameters
[nsets npar] = size(controls.large);

%% initialize some variables
controls.keepers = true(nsets,npar); % initialize flag matrix with ones
controls.small = [];
bi = isnan(controls.bias);

%% estimate B instead of z (experimental)
% Set this variable to false for normal behavior. Setting it to true causes
% DMAT to estimate B = z/a instead of z. The fourth element in the
% DesignMatrix field will then be applied to B, not z. Don't forget to
% switch it back when you're done with it! Also, you need to make the same
% change in BIGGER.M.
b_not_z = controls.fitbnotz;

% Possibilities for future implementation:
% b_not_z = getpref('dmatoolbox','version')==0.41; % For EJ
% b_not_z = getpref('dmatoolbox','estimating')=='b';
% b_not_z = controls.estimate_b_not_z;

%% If estimating B, replace z with B
if b_not_z
    controls.large(:,4) = controls.large(:,4)./controls.large(:,1);
end
    
%% evaluate restrictions one parameter at a time
for ctr = 1:npar
    if isempty(controls.designmat{ctr}) 
        controls.designmat{ctr} = eye(nsets);
    end
    dm = controls.designmat{ctr};
    fi = isnan(controls.fixvals(:,ctr));
    
%% (1) if a design was implemented
    if ~iseye(dm)
        pm = linreg(controls.large(:,ctr),dm);
        controls.small = [controls.small;pm(:)];
        controls.keepers(length(pm)+1:nsets,ctr)=false;
        controls.fixvals(:,ctr)=NaN;
        if ctr==4, controls.bias(:)=NaN; end
        
%% (2) else, if fixes were set
    elseif ~all(fi) 
        pm = controls.large(fi,ctr);
        controls.small = [controls.small;pm(:)];
        controls.keepers(~fi,ctr)=false;
        if ctr==4, controls.bias(:)=NaN; end
        
%% (3) else, if biases were requested
    elseif ctr==4 && ~all(bi) 
        pm = controls.large(bi,4);
        controls.small = [controls.small;pm(:)];
        controls.keepers(~bi,ctr)=false;

%% (4) else, no restrictions
    else 
        pm = controls.large(:,ctr);
        controls.small = [controls.small;pm(:)];
    end
end