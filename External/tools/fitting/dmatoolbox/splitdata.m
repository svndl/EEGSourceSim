function dcell = splitdata(data,splitreply)
%SPLITDATA  Splits a dataset by putting each condition into a cell
%   DATACELL = SPLITDATA(DATA), where DATA is a regular N-by-3 data matrix
%   [condition response seconds], returns DATACELL, a cell matrix with one
%   cell for each condition, and a [response seconds] data set in each.
%
%   See also PROCESSDATA.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<2
    splitreply = false;
end

if isvaliddataset(data,2)
    dcell{1} = data;
elseif isvaliddataset(data,3)
    cond = unique(data(:,1));
    ncond = length(cond);
    if splitreply
        dcell = cell(ncond,2);
        for ctr = 1:ncond
            dcell{ctr,1} = data(data(:,1)==cond(ctr)&data(:,2)==0,3);
            dcell{ctr,2} = data(data(:,1)==cond(ctr)&data(:,2)==1,3);
        end
    else
        dcell = cell(ncond,1);
        for ctr = 1:ncond
            dcell{ctr} = data(data(:,1)==cond(ctr),[2 3]);
        end
    end
else
    error('DMAT:splitdata:invalidDataSet','Data set is invalid.');
end