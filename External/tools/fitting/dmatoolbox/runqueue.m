function output = runqueue(data,queue,guessorder)
% RUNQUEUE  Runs a queue of DMA models
%   OUTPUT = RUNQUEUE(DATA,OPTIONS), where DATA is a valid N-by-3 data set
%   and OPTIONS is a R-by-Q structure array where each element is a
%   user-designed DMAT options structure, provides OUTPUT, an R-by-Q array
%   of output structures from MULTIESTV4.
%
%   OUTPUT = RUNQUEUE(DATA,OPTIONS,GUESSES), where GUESSORDER is an
%   R-by-1 vector of integers that index the OPTIONS structure, initiates
%   the estimation of OPTIONS(C) with the optimum of OPTIONS(GUESSES(C)).
%   Default is GUESSES = 0:(R*Q)-1.
%
%   See also MULTIESTV4, Q2M.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Now prints model name while estimating (useful for long queues
%  running on remote computers).

if ~isvaliddataset(data)
    error('DMAT:runqueue:badDataSet','Input is not a valid data set.')
end

[r c]=size(queue);
n = r*c;
if nargin<3,guessorder = 0:n-1; end

trysaving = 1;
try
    save('~$lastqueue','queue')
catch
    warning('DMAT:runqueue:noWritingPrivileges',...
        ['Runtime back-up procedure failed.\nYou are advised to check ',...
        'your disk access rights.\nIt is safer to run DMAT in a direct',...
        'ory where it can make runtime back-ups.'])
    trysaving = 0;
end
for ctr=1:n
    if ctr>1 && (~isfield(queue(ctr),'Guess') || ...
            isempty(queue(ctr).Guess))
        queue(ctr).Guess = output(guessorder(ctr)).Minimum;
    end
    if isfield(queue(ctr),'Name')
        fprintf('    Starting model %i of %i: %s\n',ctr,n,queue(ctr).Name);
    end
    output(ctr) = multiestv4(data,queue(ctr));
    if trysaving,save('~$queuetemp','output','data'),end
end

output = reshape(output,r,c);