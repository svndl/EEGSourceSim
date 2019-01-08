function varargout = dmatremove(varargin)
%DMATREMOVE  Removes the DMA Toolbox
%   DMATREMOVE() permanently removes all files and settings of the DMA
%   Toolbox, including all files in the DMATOOLBOX folder.
%
%   DMATREMOVE('I really do') does not ask for confirmation or print
%   warnings.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Prepared for possible (future) interaction with INSTALLER.

if ~nargin
    % ?
    fprintf('\n\nWARNING!\n')
    r = sprintf(['This will permanently remove all DMAT files, directories,',...
        'subdirectories and settings! Any files you have saved in DMAT '...
        'directories will be lost.']);
    disp(parsetoline(r))
    r = input('Type ''I really do'' if you really want to uninstall DMAT: ','s');
else
    r = varargin{1};
end

if ~strcmp(r,'I really do')
% :)
    disp 'DMAT was not uninstalled.'
    if nargout
        varargout{1} = 1;
    end
    return
end

% :'(
dms = dmatsite;
dmatdir = getpref('dmatoolbox','dmatdir');

m1 = parsetoline(sprintf(['The automated uninstaller failed to remove ',...
    'DMAT. To manually remove it, first type ''''rmpref(''dmatoolbox''',...
    ')'''' and then delete all the files in the DMAT directory: %s.\n'],...
    dmatdir));
m2 = parsetoline(sprintf(['The DMA Toolbox was uninstalled. To get the',...
    ' latest version of DMAT, please visit %s.'],dms));

try
    clear function cdfdif
    st = warning('off','MATLAB:RMDIR:RemovedFromPath');
    rmdir(dmatdir,'s')
    warning(st)
    rmpref('dmatoolbox')
    report = 1;
catch
    report = 0;
end

if ~report && ~nargin
    disp(m1)
end

% cu
if report && ~nargin
    disp ' '
    disp(m2)
    fprintf('\nKind regards,\n              - The DMAT Team\n')
end

if nargin && nargout
    varargout{1} = dmatdir;
end