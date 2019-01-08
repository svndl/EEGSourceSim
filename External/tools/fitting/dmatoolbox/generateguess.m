function options = generateguess(options,dcel)
%GENERATEGUESS  Generates an initial guess for the optimization algorithm
%   OPTIONS = GENERATEGUESS(OPTIONS,DCEL), where OPTIONS is a valid options
%   structure for DMAT, and DCEL is a split data set as returned by
%   SPLITDATA, returns OPTIONS, with the .GUESS field updated.
%
%   See also EZDIFF and PERTURB.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.3: Replaced REGRESS with a DMAT gateway LINREG to remove
%  dependence on Statistics Toolbox.
%  Edit 0.3: Added catch to avoid multiple repeats of
%  DMAT:GenerateGuess:BadGuess.

%% First, if no guess was given, generate one from EZDIFF
if isempty(options.Guess)
    % If necessary, generate a guess based on the EZDIFF algorithm
    for ctr = 1:options.nsets
        [v,a,Ter] = ezdiff(dcel{ctr});
        if Ter<=0, Ter = .3; end % sometimes ezdiff gives Ter<0
        options.Guess(ctr,:) = [a Ter .2 a/2 .45*a .9*Ter v];
        % eta, sZ & st are arbitrary guesses
    end
    if options.GuessMethodScalar == 2
        options.Guess = perturb(options.Guess);
    end
end

%% If a mixed model is used, add guesses for pi and lambda
if options.npar==9
    if size(options.Guess,2)<8
        options.Guess(:,8:9) = .97;
    elseif size(options.Guess,2)<9
        options.Guess(:,9) = .97;
    end
else
    options.Guess=options.Guess(:,1:7);
end

%% Now, make sure guesses comply with design
for ctr=1:options.npar
    if ~isempty(options.DesignMatrix{ctr})
        pm = linreg(options.Guess(:,ctr),options.DesignMatrix{ctr});
        options.Guess(:,ctr) = options.DesignMatrix{ctr}*pm;
    end
end

%% Check if this is in the parameter space (!)
[b l]=isgood(options.Guess);
ctr=0;
while ~b
    ctr=ctr+1;
    if ctr>50
        options.Guess = repmat([.003 .001 0 .0015 0 0 0],options.nsets,1);
        if options.npar==9
            options.Guess(:,[8 9]) = .97;
        end
        [b l]=isgood(options.Guess);
        if ~b
            warning('DMAT:GenerateGuess:CantFindOne',...
                'Automatic guess generator could not locate any point in the parameter space.')
            disp('If DMAT does not find a solution, you may need to input an initial guess manually.')
            disp('The last attempted point was')
            disp(options.Guess)
            disp('Which is illegal in the following fields:')
            disp(l)
        end
        break
    end
    options.Guess(l)=options.Guess(l)/2;
    [b l]=isgood(options.Guess);
    wst = 'Automatically generated guess was outside parameter space. Generating new guess.';
    lwa = lastwarn;
    if ~strcmp(lwa,wst)
        warning('DMAT:GenerateGuess:BadGuess',wst)
    end
end