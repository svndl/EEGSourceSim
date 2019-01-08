function out = eachcell(fcn,in,unwise)
% EACHCELL  Performs a function on each cell of a cell matrix.
%  OUT = EACHCELL(FCN,IN), performs function FCN (a function handle or a
%  quoted string containing the name of a function) on each cell of cell
%  matrix IN, and puts the result in the corresponding cell of OUT.
%  EACHCELL gives an error if the function FCN fails on any element.
%
%  OUT = EACHCELL(FCN,IN,UNWISE), where UNWISE is 1 (or true) does not
%  return such an error, but places a NaN in the corresponding cell.
%
%  See also CELLFUN.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<3
    if nargin<2
        error('DMAT:eachcell:NotEnoughInputs',...
            'EACHCELL requires at least 2 input arguments.')
    end
    unwise = false;
end

if ~iscell(in)
    error('DMAT:eachcell:BadInput',...
        'The second input argument for EACHCELL must be a cell matrix.')
end

sz = size(in);
out = cell(sz);

for a=1:prod(sz)
    try
        out{a} = feval(fcn,in{a});
    catch
        if unwise
            out{a} = NaN;
        else
            error('DMAT:eachcell:FunctionDidntWork',...
                ['Input function did not work on element with linear ',...
                'index %i.\nError was:\n\n   %s'],a,lasterr)
        end
    end
end