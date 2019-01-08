function varargout = namepars(varargin)
%NAMEPARS  Shows and names free parameters
%   C = NAMEPARS(OUTPUT), where OUTPUT is a DMAT output structure, returns
%   C, a k-by-3 cell matrix containing parameters' names in the first
%   column, their index in the design parameter matrix in the second
%   column, and their estimate in the third column.
%
%   NAMEPARS(OUTPUT) simply prints a list to the screen.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if ~nargin
    error('DMAT:namepars:notEnoughInputs',...
        'NAMEPARS requires at least one input variable.')
end

output = varargin{1};
if ~isstruct(output)
    error('DMAT:namepars:badOutputStructure',...
        'OUTPUT should be a valid output structure.')
end

if numel(output)>1
    for ctr=1:numel(output)
        if nargout
            varargout{1}{ctr} = namepars(output(ctr));
        else
            fprintf('Model ''%s'':\n',output(ctr).Name)
            namepars(output(ctr),varargin{2:end})
        end
    end
    return
end

c=find(output.Options.controls.keepers);
npars = length(c);
pna = cell(npars,3);
ps = cell(npars,1);

for a=1:npars,
    [nu pa] = ind2sub(size(output.Options.controls.keepers),c(a));
    switch pa
        case 1, pn = 'a';
        case 2, pn = 'Ter';
        case 3, pn = 'eta';
        case 4, pn = 'z';
        case 5, pn = 'sZ';
        case 6, pn = 'st';
        case 7, pn = 'v';
        case 8, pn = 'pi';
        case 9, pn = 'gamma';
        otherwise, warning('DMAT:namepars:unknownParameter',...
                'There are more than 9 columns in the parameter matrix.')
    end
    pna{a,1} = pn;
    pna{a,2} = nu;
    pna{a,3} = output.Options.controls.small(a);
    ps{a} = sprintf('%s%i',pn,nu);
end

if ~nargout
    fprintf('  Parameter list:\n')
    for a=1:npars
        fprintf('%5u:%8s%10.6f\n',a,...
            [pna{a,1} '(' num2str(pna{a,2}) ')'],pna{a,3})
    end
else
    if nargin<2 || ~varargin{2}
        varargout{1} = pna;
    else
        varargout{1} = ps;
    end
end
