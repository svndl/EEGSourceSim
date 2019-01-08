function varargout = waldtable(output)
% WALDTABLE  Produces a table with parameter estimates and Wald tests
%   WALDTABLE produces a table with the point estimate of each basic
%   parameter, its standard error of estimation, the Wald statistic Z, and
%   the p-value under a correct reference distribution, which is the
%   standard normal in the regular case. In the case where the test value
%   (given in the last column) is a boundary value, the correct test uses
%   Z^2 ~ .5*Chi2(0) + .5*Chi2(1). This is the case for eta, sz, st, pi,
%   and gamma.
%
%   Usage: 
%     WALDTABLE(OUTPUT), where OUTPUT is an output structure (or an array
%     of such), displays the table.
%
%     T = WALDTABLE(OUTPUT), returns the table to a cell T (or a nested
%     cell if OUTPUT is an array).
%
%   See also QTABLE and MODELFITTABLE.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Added in version 0.4.

if numel(output)>1
    for ctr=1:numel(output)
        if nargout
            varargout{1}{ctr} = waldtable(output(ctr));
        else
            waldtable(output(ctr))
            disp ' '
        end
    end
    return
end

if isfield(output,'Name')
    fprintf('   Wald table for model ''%s''\n',output.Name)
end
pns = namepars(output);
bounds = cellfun(@(x)any('espg'==x(1)),pns(:,1))+1;
h0s = cellfun(@(x)any('pg'==x(1)),pns(:,1));
sees = sqrt(diag(inv(output.Hessian)));
sees(imag(sees)<=1e-8)=real(sees); % Remove small imaginary parts
sees(imag(sees)~=0)=NaN;
walds = ([pns{:,3}]'-h0s)./sees;
ps = 2*(1-normcdf(abs(walds)));
ps = ps./bounds;
pns = [pns num2cell([sees(:) walds(:) ps(:) h0s(:)])];

if ~nargout
    a=pns';
    fprintf('%9s','Par','Est','SEE','Wald','p','H0',sprintf('\n'))
    fprintf('%6s(%i)%9.4f%9.4f%9.4f%9.4f%9.4f\n',a{:})
else
    varargout{1} = pns;
end