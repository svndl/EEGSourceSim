function [stderr warn]= processhessian(hessian,options)
%PROCESSHESSIAN  Process the Hessian matrix to get standard errors
%   [STDERR WARN]= PROCESSHESSIAN(HESSIAN,OPTIONS), where OPTIONS is a
%   valid options structure for DMAT (as it is found in the output.Options
%   field!) and HESSIAN is the Hessian matrix at the minimum of the
%   deviation function, as returned by GENALG, returns STDERR, a matrix
%   with a standard error for each diffusion parameter and WARN, a cell
%   matrix with possible warning messages. 
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Moved tolerance for positive-definiteness to 1e-6 from 1e-10.

warn = {};
[nsets npar] = size(options.controls.large);

%% See if eigenvalues aren't too small
egv = eig(hessian);
if any(egv<0&egv>-1e-6)
    warning('DMAT:processhessian:smallEigenvalues',...
        'Some Hessian eigenvalues are negative but close to zero.')
    hessian = hessian + eye(size(hessian,1)).*1e-5;
    % TODO % Perhaps the constant added could be larger?
end

%% check Hessian properties
if rank(hessian)<npar
    wstr='Hessian is not of full rank.';
    warn{end+1} = wstr;
    warning('DMAT:processhessian:HessianNotOfFullRank',wstr)
end
if ~isposdef(hessian,1e-6)
    wstr = 'Hessian is not positive definite.';
    warn{end+1} = wstr;
    warning('DMAT:processhessian:HessianNotPositiveDefinite',wstr)
end

%% prepare some variables
covarmat = inv(hessian); % covariance matrix
normvar = diag(covarmat); % variances if no design
varia = zeros(nsets,npar); % allocate for variances
matctr = 1; % keep track of free parameters

%% evaluate restrictions one parameter at a time
for ctr = 1:npar
    dm = options.controls.designmat{ctr};
    fi = ~isnan(options.controls.fixvals(:,ctr));
    bi = ~isnan(options.controls.bias);
    
%% (1) if a design was implemented
    if ~iseye(dm)   
        nfree = size(dm,2);
        indcov = matctr:matctr+nfree-1;
        weird = dm*covarmat(indcov,indcov)*dm';
        varia(:,ctr) = diag(weird);
        matctr = matctr+nfree;

%% (2) else, if fixes were set
    elseif any(fi) 
        nfree = sum(~fi);
        varia(~fi,ctr) = normvar(matctr:matctr+nfree-1);
        varia( fi,ctr) = 0;    
        matctr = matctr+nfree;
                
%% (3) else, if biases were requested
    elseif ctr==4 && any(bi) 
        nfree = sum(~bi);
        varia(~bi,4) = normvar(matctr:matctr+nfree-1);
        varia( bi,4) = varia( bi,1).*(options.controls.bias(bi).^2);
        matctr = matctr+nfree; 
                
%% (4) else, no restrictions
    else 
        nfree = nsets;
        varia(:,ctr) = normvar(matctr:matctr+nfree-1);
        matctr = matctr+nfree;
    end
end
%%

stderr = sqrt(varia);