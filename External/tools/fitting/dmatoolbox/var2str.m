function str = var2str(var)
% VAR2STR  Print variable contents to a MATLAB string
%   STR = VAR2STR(VAR), where VAR is any standard matlab variable, returns
%   STR, a character string that MATLAB would read as VAR. (It performs
%   essentially the opposite function of EVAL.)
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

w = whos('var');
switch w.class
    case {'logical','double'}
        str = m2s(var);
    case {'int8','uint8','int16','uint16',...
            'int32','uint32','int64','uint64','single'}
        str = sprintf('%s(%s)',w.class,m2s(var));
    case 'cell'
        str = c2s(var);
    case 'struct'
        str = struct2string(var);
    case 'function_handle'
        str = func2str(var);
    case 'char'
        str = ['''' var ''''];
    otherwise
        error('DMAT:var2str:unknownClass',...
            'Unknown variable type ''%s''.', w.class)
end
        
function s = c2s(m)
% C2S  Print cell contents to a MATLAB string

s='{'; [r c]=size(m);
for ri = 1:r
    for ci = 1:c-1
        s=sprintf('%s%s,',s,var2str(m{ri,ci}));
    end
    if ri<r
        s=sprintf('%s%s;',s,var2str(m{ri,c}));
    else
        s=sprintf('%s%s}',s,var2str(m{ri,c}));
    end
end


function s=m2s(m)
% M2S  Print double array contents to a MATLAB string

[r c]=size(m);

% Catch empty matrix
if r<1 || c<1
    s = '[]';
    return
end

% Catch strings
if ischar(m)
    s=sprintf('''%s''',m);
    return
end

% Catch singleton
if r*c==1
    s = sprintf('%g',m);
    return
end

% Otherwise, grow string
s='[';
for ri = 1:r
    for ci = 1:c-1
        s=sprintf('%s%g,',s,m(ri,ci));
    end
    if ri<r
        s=sprintf('%s%g;',s,m(ri,c));
    else
        s=sprintf('%s%g]',s,m(ri,c));
    end
end


function string = struct2string(structure)
% STRUCT2STRING  Print array contents to a MATLAB string

fnm = fieldnames(structure);
string = 'struct(';
l = length(fnm);
for a = 1:l
    if a==1
        string = [string '''' fnm{a} ''',' ...
            var2str(structure.(fnm{a})) ',...\n'];
    elseif a<l
        string = [string '   ''' fnm{a} ''',' ...
            var2str(structure.(fnm{a})) ',...\n'];
    else
        string = sprintf([string '   ''' fnm{l} ''',' ...
            var2str(structure.(fnm{l})) ')']);
    end
end