function varargout = paramplot(param,errbar,handle)
%PARAMPLOT  Plot parameters over conditions
%   PARAMPLOT(PARAM,ERRBAR,HANDLE) plots the values in PARAM with error
%   bars determined by the values in ERRBAR. These two should be vectors of
%   the same size. The handle in HANDLE should be to an AXES object.
%
%   PARAMPLOT(PARAM,ERRBAR) makes a new AXES object.
%
%   H = PARAMPLOT(PARAM,ERRBAR) returns a handle to the AXES object used.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<1
    error('DMAT:paramplot:notEnoughInputs',...
        'PARAMPLOT requires at least one input.')
end

[pr pc] = size(param);
if pr~=1 && pc~=1
    error('DMAT:paramplot:parameterMustBePsiVector',...
        'PARAM should be a vector.')
end

Marker = 'none';
MarkerEdgeColor = 'auto';

if nargin<2
    errbar = zeros(size(param));
elseif ~all([pr pc]==size(errbar))
    error('DMAT:paramplot:errbarSizeNotRight',...
        'ERRBAR should be the same size as PARAM.')
elseif any(~isreal(errbar)) || any(isinf(errbar)) || any(isnan(errbar))
    errbar(:)=0;
    Marker = 'o';
    MarkerEdgeColor = [1 0 0];
end

if nargin<3 || ~ishandle(handle) || ~strcmp(get(handle,'type'),'axes')
    if nargin>2 % handle was provided but not good
        warning('DMAT:paramplot:badHandleForAxes',...
            'Third argument was not a handle to an AXES object.')
    end
    handle = axes;
end

cond = 1:length(param);
h = errorbar(cond,param,errbar);
set(h,'Parent',handle,'Marker',Marker,'MarkerEdgeColor',MarkerEdgeColor);

if nargout
    varargout{1} = handle;
end