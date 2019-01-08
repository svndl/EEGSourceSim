function varargout = fitlast(par,output)
%FITLAST  Calculate the fit of any given parameter set to the last model
%   D = FITLAST(PAR,OUTPUT), where PAR is a full parameter set, and OUTPUT
%   is the output from a DMAT run, returns D, the deviance measure for that
%   parameter set.
%
%   FITLAST(PAR,OUTPUT), prints the fit to the screen.
%   
%   See also MULTIESTV4, PLOTPARREG.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

output(1).Options.controls.large = par;
output(1).Options.controls = smaller(output(1).Options.controls);
D = output(1).Options.objecfun(output(1).Options.controls.small);

if nargout
    varargout{1} = D;
else
    if ismember(output(1).Options.EstimationMethodScalar,2:4)
        meth = 'X2';
    elseif ismember(output(1).Options.EstimationMethodScalar,5:7)
        meth = 'MLF';
    else
        meth = 'Deviance';
    end
    fprintf(1,'\nBadness-of-this: %s = %f\n\n',meth,D);
end