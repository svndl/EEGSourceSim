function varargout=isconsistent(controls)
%ISCONSISTENT  Check if model restrictions do not overlap
%   ISCONSISTENT(CONTROLS), where CONTROLS is a DMAT controls structure,
%   gives warnings if any one parameter has more than one type of
%   restriction (design, fixed parameters, or specificbiases).
%
%   B = ISCONSISTENT(CONTROLS) also returns a boolean indicating if the
%   input is consistent (true) or not (false).
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Corrected warning tag.

npars = size(controls.fixvals,2);
pars = {'a','Ter','eta','z','sz','st','v','pi','gamma'};

wtag = 'DMAT:isconsistent:overfixing';
ok = true;

for ctr=1:npars
    if hasdesign
        if hasfixes && hasbiases
            warning(wtag,['Parameter ''%s'' has design restrictions, '...
                'but also has fixed parameters and specific biases. '...
                'Parameter fixes and specific biases will be ignored.'],...
                pars{ctr})
            ok = false;
        elseif hasfixes
            warning(wtag,['Parameter ''%s'' has design restrictions, '...
                'but also has fixed parameters. Parameter fixes will '...
                'be ignored.'],pars{ctr})
            ok = false;
        elseif hasbiases
            warning(wtag,['Parameter ''%s'' has design restrictions, '...
                'but also has specific biases. Specific biases will '...
                'be ignored.'],pars{ctr})
            ok = false;
        end
    elseif hasfixes && hasbiases
        warning(wtag,['Parameter ''%s'' has parameter fixes, but '...
            'also has specific biases. Specific biases will be '...
            'ignored.'],pars{ctr})
            ok = false;
    end
end

if nargout
    varargout = {ok};
end

    function c=hasdesign
        c = ~isempty(controls.designmat) && ...
            ~isempty(controls.designmat{ctr}) && ...
            ~iseye(controls.designmat{ctr});
    end

    function c=hasfixes
        c = ~isempty(controls.fixvals) && ...
            ~isempty(controls.fixvals(:,ctr)) && ...
            ~all(isnan(controls.fixvals(:,ctr)));
    end

    function c=hasbiases        
        c = ctr==4 && ...
            ~isempty(controls.bias) && ...
            ~all(isnan(controls.bias));
    end

end