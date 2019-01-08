function output=m2q(fname,silent_flag)
% M2Q  Read an M-file produced by DMATGUI and extract OPTIONS
%   OPTIONS = M2Q(FILENAME), where filename is the full path to an M-file
%   which was produced by the 'save to M-file' function in DMATGUI,
%   extracts the OPTIONS array (the queue) from that M-file.
%   The reading and interpreting of an M-file is a bit tricky, so errors
%   might occur. In general, if the input file is an *unedited* M-file from
%   DMATGUI, it should work. If the M-file has been significantly altered,
%   recovering the OPTIONS structure might not work.
%
%   In practice, this function reads and executes the script until the
%   function 'runqueue' is encountered (then it stops), and then tries to
%   find a structure that is called 'options' or has the field 'Name' in
%   it. Errors will occur if no structures are found, more than one is
%   found but none are called 'options' or have a 'Name' field. If more
%   than one with that field is found, the first one encountered will be
%   chosen.
%
%   The function will not work if there are lines of code (except comments)
%   in the M-file that do not end with ';', ',', or '...'.
%
%   See also DMATGUI.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Removed dependence on KEEP.

if nargin<2
    if nargin<1
        error('DMAT:m2q:notEnoughInputs',...
            'M2Q requires at least one input.')
    end
    silent_flag = false;
end

m2q_structure.input.silent = silent_flag;
m2q_structure.input.fname = fname;
m2q_structure.output.msg = [];

if ~strcmp(m2q_structure.input.fname(end-1:end),'.m')
    m2q_structure.input.fname = [m2q_structure.input.fname '.m'];
end

fid = fopen(m2q_structure.input.fname,'r');
mtext = textscan(fid,'%s','delimiter','\n');
fclose(fid);

mtext = mtext{1};
nlines = length(mtext);
commstr = [];

for l=1:nlines
    str = char(mtext{l});
    if length(str)>0&&str(1)~='%'
        if isempty(strfind(str, 'runqueue'))
            commstr = [commstr str];
        else
            break
        end
    end
end

l = findstr(commstr,'...');
l = [l l+1 l+2];
commstr(l)=[];

clear fid fname l mtext nlines silent_flag str

try
    evalc(commstr);
    if exist('fn','var')
        m2q_structure.output.filename = fn;
    else
        m2q_structure.output.filename = [];
    end
catch
    output = struct('msg',sprintf('Couldn''t read queue from %s.',...
        m2q_structure.input.fname),...
        'filename',[],...
        'queue',[]);
    if ~m2q_structure.input.silent
        error('DMAT:m2q:CouldntReadQueue',m2q_structure.output.msg);
    end
    return
end

vars = whos;
structures = ~cellfun('isempty',strfind({vars.class},'struct'));
if sum(structures)>1
    nms={vars.name};
    iscalledoptions = cell2mat(eachcell(@(x) strcmp(x,'options'),nms));
    if any(iscalledoptions)
        choice = find(iscalledoptions);
    else
        choice = 0;
        for ctr = 1:length(vars)
            if structures(ctr) && ...
                    eval(sprintf('isfield(%s,''Name'');',nms{ctr}));
                choice=ctr;
                break
            end
        end
    end
    if choice
        evalc('options = %s;',vars(choice).name);
    else
        m2q_structure.output.msg = ...
            'Multiple structures were found in the M-file, can''t choose.';
        if ~m2q_structure.input.silent
            error('DMAT:m2q:cantChoose',m2q_structure.output.msg)
        end
    end
elseif sum(structures)<1
    m2q_structure.output.msg = 'No structures found in script.';
    if ~m2q_structure.input.silent
        error('DMAT:m2q:noStructures',m2q_structure.output.msg)
    end
else
    evalc('options = %s;',vars(structures).name);
end

if ~isempty(m2q_structure.output.msg)
    m2q_structure.output.queue = [];
    m2q_structure.output.filename = [];
else
    output = struct('msg',m2q_structure.output.msg,...
        'filename',m2q_structure.output.filename,...
        'queue',options);
end