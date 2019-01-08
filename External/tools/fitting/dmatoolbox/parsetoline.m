function ml=parsetoline(sl,maxlen,fl)
%PARSETOLINE  Parses a string so it fits on a screen
%   MULTILINE = PARSETOLINE(SINGLELINE,MAXLENGTH), where SINGLELINE is a
%   text string and MAXLEN is an integer indicating the maximum horizontal
%   lenght of a string, returns MULTILINE, a string that is parsed to fit
%   on such a line. Parsing takes place only at spaces.
%
%   MULTILINE = PARSETOLINE(SINGLELINE,MAXLENGTH,FILL), pads the string
%   with spaces so that each line contains exactly MAXLENGTH characters.
%   (If you do not want to provide MAXLENGTH, use [] as a placeholder.)
%
%   MULTILINE = PARSETOLINE(SINGLELINE) parses to fit to the current
%   command window size.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if nargin<3
    fl = 0;
    if nargin<2 || isempty(maxlen)
        cwsz = get(0,'CommandWindowSize');
        maxlen = cwsz(1);
    end
end

if iscell(sl)
    ml = cellfun(@(x)parsetoline(x,maxlen,fl),sl,'UniformOutput',0);
    return
end

if ~ischar(sl)
    error('DMAT:parsetoline:notAString',...
        'First input to PARSETOLINE should be a character array.')
end

if length(sl)<=maxlen
    ml = sl;
    if fl
        ml = [sl repmat(' ',1,maxlen-length(sl))];
    end
    return
end

ind = findstr(' ',sl(1:maxlen));
if isempty(ind)
    warning('DMAT:parsetoline:noSpacesFound',...
        ['Not enough spaces found to parse string, longest word is ',...
        '%i characters.'], max(diff([0 find(isspace(sl))]))-1)
    ml = sl;
else
    fp = sl(1:ind(end));
    if fl
        fp = [fp repmat(' ',1,maxlen-length(fp))];
    end
    ml = sprintf('%s\n%s',fp,parsetoline(sl(ind(end)+1:end),maxlen,fl));
end
