function varargout=loaddmatdata(fn)
%LOADDMATDATA  Load a data file for DMAT
%   DATA = LOADDMATDATA(FILENAME) loads the data found in FILENAME into the
%   matrix DATA. This function works with these formats:
%       Extension            Properties                    Returns
%            .txt       Space delimited              Full contents
%            .csv       Comma delimited              Full contents
%       .tab .dat         Tab delimited              Full contents
%            .mat           MATLAB file       Contents of variable
%                                               'data' or first in
%                                               alphabetical order
%   Regardless of file type, LOADDMATDATA will return an error if the
%   loaded data set is not in a proper format for DMAT.
%
%   [DATA ERRORMSG WARNINGMSG] = LOADDMATDATA(FILENAME) does not throw
%   errors or warnings to the base workspace, but passes them as output to
%   the caller. ERRORMSG and WARNINGMSG are then both 2-by-1 cell arrays,
%   the first element containing an identifier string and the second a text
%   string. If there were no errors or warnings, they are empty.
%
%   LOADDMATDATA sorts the loaded datafile in ascending order of condition,
%   descending of response, and ascending of RT (in that order).
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%% Check call format
if ~nargin
    error('DMAT:loaddmatdata:notEnoughInputs',...
        'LOADDMATDATA requires at least one input variable.')
end

if nargout~=1 && nargout~=3
    error('DMAT:loaddmatdata:inappropriateNumberOfOutputs',...
        'LOADDMATDATA gives one or three variables as output.')
end

%% Initialise error/warn flags
err = 0;
data = [];
errormsg = [];
errortag = [];
warnmsg = [];
warntag = [];

%% Load data from file
ext = fn(end-2:end);
switch ext
    case {'txt'}
        try
            data = dlmread(fn,' ');
        catch
            err = 1;
            [errormsg errortag] = lasterr;
        end
    case {'tab';'dat'}
        try
            data = dlmread(fn,'\t');
        catch
            err = 1;
            [errormsg errortag] = lasterr;
        end
    case {'mat'}
        try
            evalc('s = load(fn);');
            if isfield(s,'data')
                data = s.data;
            else
                name = sort(fieldnames(s));
                data = s.(name{1});
            end
        catch
            err = 1;
            [errormsg errortag] = lasterr;
        end
    case {'csv'}
        try
            data = dlmread(fn,',');
        catch
            err = 1;
            [errormsg errortag] = lasterr;
        end
    otherwise
        err = 1;
        errortag = 'DMAT:loaddmatdata:unknownFileType';
        errormsg = 'Unknown file type';
end

%% If loading successful, check if data set is valid
if ~err
    n = size(data,1);
    [isd3 erd3] = isvaliddataset(data,3);
    valid = isd3;
    if ~isd3 % if it's a good 3-col data set, skip along, otherwise...
        [isd2 erd2] = isvaliddataset(data,2);
        valid = isd2;
        if isd2 % if it's a good 2-col data set add first col and go on
            data = [ones(n,1) data];
        elseif ~isd3 && ~isd2 % otherwise give up, report what went wrong
            data = [];
            c = size(data,2);
            if c==3
                errormsg = erd3;
            elseif c==2
                errormsg = erd2;
            else
                errormsg = ['Data set does not have a right number of '...
                    'columns.'];
            end
            errormsg = ['Loading failed: ' errormsg];
            errortag = 'DMAT:loaddmatdata:invalidDataSet';
        end
    end
end

%% Check if suitable for DMA
if ~err && valid
    md = mean(data(:,3));
    if md>10
        warntag = 'DMAT:dmatgui:tooHighRTs';
        warnmsg = ['The mean response time is ' num2str(md,'%4.3f'),...
            ' seconds, which is too high for diffusion model analysis.'];
    end
end

%% Prepare output
% if three output arguments -> pass errors and warnings to caller
% otherwise -> throw to base
if nargout == 1
    if isempty(errormsg)
        varargout{1} = sortrows(data,[1 -2 3]);
    end
    if ~isempty(errormsg)
        error(errortag,errormsg)
    elseif ~isempty(warnmsg)
        warning(warntag,warnmsg)
    end
elseif nargout == 3
    if ~isempty(errormsg)
        varargout{2} = {errortag;errormsg};
    else
        varargout{1} = sortrows(data,[1 -2 3]);
        varargout{2} = [];
    end
    if ~isempty(warnmsg)
        varargout{3} = {warntag;warnmsg};
    else
        varargout{3} = [];
    end
end

