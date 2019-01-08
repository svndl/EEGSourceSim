function [a bad] = isgood(par)
%ISGOOD  Check if parameter set is inside the general parameter space
%   [GOOD BAD] = ISGOOD(PAR), where PAR is any parameter set, returns GOOD,
%   a logical scalar which is true if the parameter set is good, and false
%   otherwise. It also returns BAD, which is a logical matrix indicating
%   for each parameter whether or not it is good.
%
%   A parameter is 'good' if it is in its possible range and 'bad'
%   otherwise. A tolerance of .001 is built in. Drift rates have maxima of
%   5 and minima of -5, boundary separations have a maximum of 2, and etas
%   have maxima of 0.5.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.3: Now uses REPMAT for speed.
%   Edit 0.4: a now limited to 0 < a < 2

[ns np]=size(par);
if np~=7&&np~=9
    a = false;
    bad = ones(ns,np);
    return
end

maxs = [repmat(2,ns,1), ...                        % max for a
        repmat(inf,ns,1), ...                      % max for Ter
        repmat(.5,ns,1), ...                       % max for eta
        par(:,1)-.001, ...                         % max for z0
        2*min(par(:,1)-par(:,4),par(:,4))-.001, ...% max for sZ
        2*par(:,2)-.001, ...                       % max for st
        repmat(5,ns,1)];                           % max for nu

mins = [par(:,4)+.001,...                 % min for a
        repmat(0,ns,5), ...               % min for Ter, eta, z0, sZ, st
        repmat(-5,ns,1)];                 % min for nu

bad = (par(:,1:7)<mins | par(:,1:7)>maxs | isnan(par(:,1:7)));
a = any(bad(:));

if np==9
    bad = (par(:,8:9)<0 | par(:,8:9)>1);
    a = a|any(bad(:));
end

a=~a;

% if ~a,keyboard,end