function output = runqueue_q(data,queue,guessorder)
% RUNQUEUE_Q  Runs a queue of DMA models with quantile input
%   RUNQUEUE_Q works just like RUNQUEUE, but takes the same data input as
%   QUANTEST.
%
%   See also QUANTEST.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Now catches saving error if user has no write access.
%  Edit 0.4: Now prints model name while estimating (useful for long queues
%  running on remote computers).

%  Forward compatibility notice: This function will likely be absorbed into
%  RUNQUEUE / MULTIESTV4 in a future DMAT release.

[r c]=size(queue);
n = r*c;
if nargin<3,guessorder = 0:n-1; end

trysaving = 1;
try
    save('~$lastqueue','queue')
catch
    warning('DMAT:runqueue_q:noWritingPrivileges',...
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
    output(ctr) = quantest(data,queue(ctr));
    if trysaving,save('~$queuetemp','output','data'),end
end

output = reshape(output,r,c);