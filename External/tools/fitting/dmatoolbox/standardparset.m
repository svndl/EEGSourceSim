function param=standardparset(parset,v)
% STANDARDPARSET  Standard diffusion model parameter sets
%   PARAM = STANDARDPARSET(X), generates standard parameter set X. You can
%   define your own if you like.
%
%   PARAM = STANDARDPARSET(0,D) (with X=zero) generates the parameter set
%                    a  Ter  eta   z  sz  st   v
%    Condition 1   .08  .30  .08 .04 .02 .02   D
%    Condition 2   .08  .30  .08 .04 .02 .02 .00
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

switch parset
    case 0
        param = [.08 .30 .08 .04 .02 .02 v
                 .08 .30 .08 .04 .02 .02 .00];
    case 1
        param = [.08 .30 .08 .04 .02 .02 .40
                 .08 .30 .08 .04 .02 .02 .25
                 .08 .30 .08 .04 .02 .02 .10
                 .08 .30 .08 .04 .02 .02 .00];
    case 2
        param = [.08 .30 .16 .04 .02 .02 .40
                 .08 .30 .16 .04 .02 .02 .25
                 .08 .30 .16 .04 .02 .02 .10
                 .08 .30 .16 .04 .02 .02 .00];
    case 3
        param = [.16 .30 .08 .08 .02 .10 .30
                 .16 .30 .08 .08 .02 .10 .20
                 .16 .30 .08 .08 .02 .10 .10
                 .16 .30 .08 .08 .02 .10 .00];
    case 4
        param = [.16 .30 .16 .08 .02 .10 .30
                 .16 .30 .16 .08 .02 .10 .20
                 .16 .30 .16 .08 .02 .10 .10
                 .16 .30 .16 .08 .02 .10 .00];
    case 5
        param = [.16 .30 .08 .08 .10 .10 .30
                 .16 .30 .08 .08 .10 .10 .20
                 .16 .30 .08 .08 .10 .10 .10
                 .16 .30 .08 .08 .10 .10 .00];
    case 6
        param = [.16 .30 .16 .08 .10 .10 .30
                 .16 .30 .16 .08 .10 .10 .20
                 .16 .30 .16 .08 .10 .10 .10
                 .16 .30 .16 .08 .10 .10 .00];
    case 7
        param = [.08 .30 .08 .03 .02 .02 -.30
                 .08 .30 .08 .03 .02 .02 -.10
                 .08 .30 .08 .03 .02 .02 .10
                 .08 .30 .08 .03 .02 .02 .30];
    case 8
        param = [.08 .30 .16 .03 .02 .02 -.30
                 .08 .30 .16 .03 .02 .02 -.10
                 .08 .30 .16 .03 .02 .02 .10
                 .08 .30 .16 .03 .02 .02 .30];
    case 9
        param = [.16 .30 .08 .06 .02 .10 -.30
                 .16 .30 .08 .06 .02 .10 -.10
                 .16 .30 .08 .06 .02 .10 .10
                 .16 .30 .08 .06 .02 .10 .30];
    case 10
        param = [.16 .30 .16 .06 .02 .10 -.30
                 .16 .30 .16 .06 .02 .10 -.10
                 .16 .30 .16 .06 .02 .10 .10
                 .16 .30 .16 .06 .02 .10 .30];
    case 11
        param = [.16 .30 .08 .06 .10 .10 -.30
                 .16 .30 .08 .06 .10 .10 -.10
                 .16 .30 .08 .06 .10 .10 .10
                 .16 .30 .08 .06 .10 .10 .30];
    case 12
        param = [.16 .30 .16 .06 .10 .10 -.30
                 .16 .30 .16 .06 .10 .10 -.10
                 .16 .30 .16 .06 .10 .10 .10
                 .16 .30 .16 .06 .10 .10 .30];
    case 13
        param = [.08 .30 .00 .04 .00 .00 .001
                 .08 .30 .00 .04 .00 .00 .10
                 .08 .30 .00 .04 .00 .00 .25
                 .08 .30 .00 .04 .00 .00 .40];
    case 14
        param = [.08 .30 .00 .04 .00 .00 .00];
    case 15
        param = [.08 .30 .16 .04 .03 .10 .10
                 .10 .30 .16 .05 .03 .10 .10
                 .12 .30 .16 .06 .03 .10 .10
                 .14 .30 .16 .07 .03 .10 .10];
    case 16
        param = [.12 .30 .16 .036 .03 .10 .10
                 .12 .30 .16 .048 .03 .10 .10
                 .12 .30 .16 .060 .03 .10 .10
                 .12 .30 .16 .072 .03 .10 .10];
    otherwise
        error('No parameter set %i is defined.',parset)
end
