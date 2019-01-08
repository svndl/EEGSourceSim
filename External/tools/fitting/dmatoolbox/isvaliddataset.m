function [b mssg] = isvaliddataset(data,cols)
% ISVALIDDATASET  Check if input is a valid data set for DMAT
%   B = ISVALIDDATASET(DATA,COLS), where DATA is the matrix to be checked,
%   and COLS is 1, 2, or 3, depending on the number of columns expected.
%      If COLS is 3, the structure should be [condition response seconds].
%      If COLS is 2, the structure should be [response seconds].
%      If COLS is 1, the structure should be [seconds].
%   B indicates if this is true.
%
%   [B,TEXT] = ISVALIDDATASET(DATA,COLS) also returns a text string
%   explaining what was wrong with the data matrix.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.
    
if nargin<2
    cols=3;
elseif ~ismember(cols,1:3)
    error('DMAT:isvaliddataset:badColumnsRequest',...
        'COLS should be 1, 2, or 3.')
end
b=true;
mssg = [];

if ~isnumeric(data)
    b=false;
    mssg = 'Input is not numeric.';
    return
end

[r c] = size(data);
if c~=cols
    b=false;
    mssg = 'Number of columns does not agree.';
    return
end

if r<2
    b=false;
    mssg = 'Not enough trials for a real data set.';
    return
end

for c=1:cols
    switch c
        case 1
            if any(data(:,end)<0)
                b=false;
                mssg = 'Negative reaction times found in last column.';
                return
            end
        case 2
            if any(~ismember(data(:,end-1),[0 1]))
                b=false;
                mssg = ['Responses other than 0 or 1 found in second ',...
                    'to last column.'];
                return
            end
        case 3
            return
    end
end
    
    