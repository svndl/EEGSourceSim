function [T,XX] = simuldiff(par,N)
%SIMULDIFF  Generates data according to a diffusion model
%   [T,X] = SIMULDIFF(PAR,N), where PAR is a valid parameter vector with
%   elements [a, Ter, eta, z, sz, st, v], and N is a scalar, returns T, a
%   1-by-N vector with reaction times (in seconds), and X, a 1-by-N vector
%   with binary responses.
%
%   DATA = SIMULDIFF(PAR,N), returns DATA, an N-by-2 data matrix with
%   binary responses in the first column and reaction times (in seconds) in
%   the second column.
%
%   The volatility parameter s is fixed to .1.
%
%   Drift rates are limited to the interval [-.5 .5] and their standard
%   deviation cannot exceed .30.
%
%   See also MULTISIMUL, SIMULDIFFWO.
%
%   This code is based on methods described in Tuerlinckx, F., Maris, E.,
%   Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods for
%   simulating the diffusion process. _Behavior Research Methods,
%   Instruments, & Computers, 33,_ 443-456.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command.
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

% Note: If you are now reading this code because you want to optimize it
% for speed, contact me and I can send you a (C) MEX-file that is an order
% of magnitude faster. It's not included in DMAT yet because it isn't very
% easy to compile (relies on GNU GSL). -jv 11/apr/2007

if nargin<2
    N = 250;
    if nargin<1
        error('DMAT:simuldiff:notEnoughInputs',...
            'SIMULDIFF requires at least one input variable.')
    end
end
if numel(par)~=7
    error('DMAT:simuldiff:badParameterSet1',...
        'The first argument to SIMULDIFF should contain exactly seven elements.')
end

a=par(1);
Ter=par(2);
eta=par(3);
z=par(4);
sZ=par(5);
st=par(6);
nu=par(7);
if nu<-.5 || nu>.5
    warning('DMAT:simuldiff:driftRateOutOfBounds',...
        'Drift rate is out of bounds, should be in [-.5 .5].')
    nu= sign(nu)*.5;
end
if eta>.3
    warning('DMAT:simuldiff:etaOutOfBounds',...
        'Standard deviation of drift rate is out of bounds, should be <.3.')
    eta=.3;
end

if eta==0, eta=1e-16; end

if ~isgood([a Ter eta z sZ st nu])
    error('DMAT:simuldiff:badParameterSet2','Not a good parameter set.')
end

T=zeros(1,N);
XX=zeros(1,N);

% parameter values
tau=.1;
% called sigma in BRMIC,2001
D=tau^2/2;

% program specifications
delta=eps;
ccc = 10001;

for j=1:N
    %if j==ccc; fprintf(1,'.'); ccc=10000+ccc; end
    r1 = randn;
    mu = nu+r1*eta;
    zz = z-sZ/2+sZ*rand;
    finish=0;
    totaltime=0;
    startpos=0;
    Aupper=a-zz;
    Alower=-zz;
    radius=min([abs(Aupper) abs(Alower)]);
    while finish==0
        lambda = 0.25*mu^2/D + 0.25*D*pi^2/radius^2;
        % eq. formula (13) in BRMIC,2001 with D = sigma^2/2 and radius = a/2
        F=D*pi/(radius*mu);
        F=F^2/(1+F^2);
        % formula p447 in BRMIC,2001
        prob=exp(radius*mu/D);
        prob=prob/(1+prob);
        dir_=2*(rand<prob)-1;
        l=-1;
        s2=0;
        while s2>l
            s2=rand; s1=rand;
            tnew=0;  told=0;
            uu=0;
            while (abs(tnew-told)>eps) || ~uu
                told=tnew;
                uu=uu+1;
                tnew = told + (2*uu+1) * (-1)^uu * s1^(F*(2*uu+1)^2);
                % infinite sum in formula (16) in BRMIC,2001
            end
            l = 1 + s1^(-F) * tnew;
            % rest of formula (16)
        end
        t = abs(log(s1))/lambda;
        % is the negative of t* in (14) in BRMIC,2001
        totaltime=totaltime+t;
        dir_=startpos+dir_*radius;
        ndrt = Ter-st/2+st*rand;
        if (dir_+delta>Aupper)
            T(j)=totaltime+ndrt;
            XX(j)=1;
            finish=1;
        elseif (dir_-delta<Alower)
            T(j)=totaltime+ndrt;
            XX(j)=0;
            finish=1;
        else
            startpos=dir_;
            radius=min(abs([Aupper Alower]-startpos));
        end
    end
end

if nargout==1;
    T = [XX' T'];
elseif ~nargout;
    fprintf(1,'\n  Simulation summary:\n\n  Mean(X) = %f\n  RT(X=1) = %f\n  RT(X=0) = %f\n\n', ...
        mean(XX), mean(T(XX==1)), mean(T(XX==0)));
else
    XX=logical(XX);
end
if ccc>10001; fprintf(1,'\n'), end