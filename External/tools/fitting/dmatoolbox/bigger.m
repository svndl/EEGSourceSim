function controls = bigger(controls)
%BIGGER  Rebuilds complete parameter set from free parameters
%   CONTROLS = BIGGER(CONTROLS)
%
%   See also SMALLER.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.31: Prepared for making B = z/a estimable in a future release.
%   Edit 0.4: Made B = z/a estimable.

%% determine conditions and parameters
[nsets npar] = size(controls.large);

%% initialize some variables
bi = isnan(controls.bias);
sm = controls.small(:);

%% estimate B instead of z (experimental)
% Set this variable to false for normal behavior. Setting it to true causes
% DMAT to estimate B = z/a instead of z. The fourth element in the
% DesignMatrix field will then be applied to B, not z. Don't forget to
% switch it back when you're done with it! Also, you need to make the same
% change in SMALLER.M.
b_not_z = controls.fitbnotz; 

% Possibilities for future implementation:
% b_not_z = getpref('dmatoolbox','version')==0.41; % For EJ
% b_not_z = getpref('dmatoolbox','estimating')=='b';
% b_not_z = controls.estimate_b_not_z;

%% evaluate restrictions one parameter at a time
for ctr = 1:npar
    dm = controls.designmat{ctr};
    fi = isnan(controls.fixvals(:,ctr));
    
%% (1) if a design was implemented
    if ~(isempty(dm)||iseye(dm)) 
        pte = size(dm,2); % number of parameters to extract
        controls.large(:,ctr) = dm*sm(1:pte); % apply
        sm(1:pte)=[]; % erase parameters
        
%% (2) else, if fixes were set
    elseif ~all(fi) 
        pte = sum(fi); % number of parameters to extract
        controls.large(fi,ctr) = sm(1:pte); % apply
        controls.large(~fi,ctr) = controls.fixvals(~fi,ctr);
        sm(1:pte)=[]; % erase parameters     
        
%% (3) else, if biases were requested
    elseif ctr==4 && ~all(bi) 
        pte = sum(bi); % number of parameters to extract
        controls.large(~bi,ctr) = controls.large(~bi,1).*controls.bias(~bi); % apply
        controls.large(bi,ctr) = sm(1:pte);
        sm(1:pte)=[]; % erase parameters
   
%% (4) else, no restrictions
    else 
        controls.large(:,ctr) = sm(1:nsets); % apply
        sm(1:nsets)=[]; % erase parameters
    end
end

%% If estimating B, replace B with z
if b_not_z
    controls.large(:,4) = controls.large(:,1).*controls.large(:,4); % replace B with z
end