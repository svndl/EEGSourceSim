function [v,a,Ter] = ezdiff(Pc,VRT,MRT)
%EZDIFF  Solves the EZDIFF system of equations
%   [V,A,TER] = EZDIFF(PC,VRT,MRT), where PC is the proportion of correct
%   responses, VRT is the variance of the correct-responses reaction time,
%   and MRT is the mean correct-responses reaction time, returns guesses
%   for V (drift rate), A (boundary separation), and TER (nondecision
%   time).
%
%   [V,A,TER] = EZDIFF(DATA), gives the same output, with DATA now an
%   N-by-2 or N-by-3 matrix with data in the format [condition response
%   seconds].
%
%   See also GENERATEGUESS, PERTURB.
%
%   The EZDIFF algorithm is due to:
%     Wagenmakers, E.-J., van der Maas, H., & Grasman, R. (in press). An
%         EZ-diffusion model for response time and accuracy. _Psychonomic
%         Bulletin & Review._
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%   Edit 0.3: Changed warning to error if non-real values pop up.
%   Edit 0.4: Changed error back to warning, simply remove imaginary part.
%   Edit 0.4: Fixed case where Pc is special but n unknown. 

vsign=1;
if nargin == 1;
    if isempty(Pc)
        VRT = 1e-9;
        Pc = .75;
        MRT = .5;
        n = 1;
    else
        dat = Pc(:,end-1:end);
        Pc = size(dat,1);
        data = dat(round(dat(:,1))==1,2);
        n = length(data);
        if n<2
            data=dat(round(dat(:,1))==0,2);
            n=length(data);
            vsign=-1;
        end
        clear dat
        VRT = var(data);
        MRT = sum(data)/n;
        Pc = n/Pc;
    end
else
    n = 100; % Good a number as any...
end

if     Pc==0.0;  Pc = 0.0 + 1/n;
elseif Pc==0.5;  Pc = 0.5 + 1/n;
elseif Pc==1.0;  Pc = 1.0 - 1/n;
end

if VRT==0.0;  VRT = 1e-9; end

s = .1;
s2 = s^2;

logit = @(x) (log(x/(1-x)));
L = logit(Pc);
x = L*(L*Pc^2-L*Pc + Pc - .5)/VRT;
v = sign(Pc-.5)*s*x^(.25); % drift rate
a = s2*logit(Pc)/v; % boundary separation
y = -v*a/s2;
MDT = (a/(2*v))*(1-exp(y))/(1+exp(y));
Ter = MRT-MDT; % nondecision time
v=v*vsign;

if ~all(isreal([v a Ter]));
    warning('DMAT:ezdiff:unlikelyParameters',...
        ['EZDIFF estimation provided unlikely parameters.\n'...
        'Pc = %f, VRT = %f, MRT = %f, yielded [v a Ter] = [%f %f %f].\n'],...
        Pc,VRT,MRT,v,a,Ter);
    v = real(v);
    Ter = real(Ter);
    a = real(a);
end