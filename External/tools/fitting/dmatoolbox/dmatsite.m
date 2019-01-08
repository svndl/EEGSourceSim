function varargout = dmatsite
%DMATSITE  The DMA Toolbox website: http://ppw.kuleuven.be/okp/dmatoolbox/
%   URL = DMATSITE() returns the web address for the DMAT website.
%
%   DMATSITE() by itself opens the website in the system browser.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

url = 'http://ppw.kuleuven.be/okp/dmatoolbox/';
if nargout
    varargout{1} = url;
else
    web(url,'-browser');
end