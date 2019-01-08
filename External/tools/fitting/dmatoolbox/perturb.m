function p = perturb(par,perc)
%PERTURB  Picks a parameter set close to another parameter set
%   NEWPAR = PERTURB(OLDPAR), returns NEWPAR, whose elements are between
%   80% and 120% of OLDPAR.
%
%   NEWPAR = PERTURB(OLDPAR,RANGE), returns NEWPAR, whose elements are
%   between OLDPAR*(1-RANGE/2) and OLDPAR*(1+RANGE/2).
%
%   See also GENERATEGUESS, EZDIFF.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<2
    if nargin<1
        error('DMAT:perturb:notEnoughInputs',...
            'PERTURB requires at least one input.')
    end
    perc=.4;
end

p = par.*(rand(size(par))*perc+(1-perc/2));
p(:,5) = .9*p(:,4);

while ~isgood(p)
    p = par.*(rand(size(par))*perc+(1-perc/2));
end