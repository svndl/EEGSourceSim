function dmatdoc(varargin)
%DMATDOC  Opens local HTML documentation for DMAT
%    DMATDOC(), without input arguments, opens the DMA Toolbox
%    documentation in MATLAB's Help Browser.
%
%    DMATDOC(LOCATION), opens it in the following programs:
%      'open' uses the MATLAB Web Browser
%      'system' uses system-dependent program (may not work on your system)
%
%    DMATDOC('web',OPTIONS), uses the MATLAB Web Browser with options:
%      '-new' opens a new MATLAB Web Browser
%      '-notoolbar' opens a MATLAB Web Browser without toolbars
%      '-noaddressbox' opens a MATLAB Web Browser without address box
%      Any combination of the above three (see WEB for help)
%      '-helpbrowser' uses the MATLAB Help Browser (default)
%      '-browser' uses your default web browser
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if ~nargin
    varargin = {'-helpbrowser'};
end
if varargin{1}(1)=='-'
    location='web';
    options = varargin;
else
    location = varargin{1};
    options = varargin(2:end);
end

url = fullfile(getpref('dmatoolbox','dmatdir'),'doc','index.html');

switch location
    case {'open','system'}
        feval(location,url);
    case 'web'
        web(url,options{:});
    otherwise
        error('DMAT:dmatdoc:unknownBrowser',...
            'LOCATION ''%s'' is unknown.',location);
end