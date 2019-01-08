function varargout = test_main(ret)
%TEST_MAIN  Diagnose basic functionalities of DMAT
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

%  Edit 0.4: Allow deployer to use this function to make new benchmark
%  file.

%% Check existance of benchmark file
tic
try
    dmsite = dmatsite;
catch
    dmsite = 'http://ppw.kuleuven.be/okp/dmatoolbox/';
end
f = filesep;
dmatdir = getpref('dmatoolbox','dmatdir');

if ~exist([dmatdir f 'test' f 'dmatdiagn.bmf'],'file')
    error('DMAT:test_main:benchmarkFileNotFound',...
        ['The file ''dmatdiagn.bmf'' was not found. Please download ',...
        'DMAT again from from the <a href="%s">DMAT website</a>.'],dmsite);
end

%% Check contents of toolbox directory
dc = dir([dmatdir f '*.*']);
foundcontents = {dc.name};
try
    dmatcontents = installer;
catch
    error('DMAT:test_main:installerNotFound',...
        ['The file ''installer.m'' was not found. Please download ',...
        'DMAT again from from the <a href="%s">DMAT website</a>.'],dmsite);
end

dmatcontents(1:2)=[];
missing = ~ismember(dmatcontents,foundcontents);
if any(missing)
    mflcell = dmatcontents(missing);
    error('DMAT:test_main:filesMissing',[parsetoline(sprintf(...
        'The files %s are missing from your DMAT installation. ',...
        [sprintf('''%s'', ',mflcell{1:end-1}) 'and ''' mflcell{end},...
        '''']),74),['You should download DMAT again from the '...
        '<a href="%s">DMAT website</a>.']],dmsite)
end

stat = warning('off','MATLAB:singularMatrix');
warning('off','DMAT:processhessian:HessianNotPositiveDefinite')
warning('off','DMAT:processhessian:HessianNotOfFullRank')

%% Diagnose each file separately
ms='DMAT diagnostics are running. This may take a few minutes.';
disp(ms)

% extract m-files
ism = cellfun(@(x)x(end)=='m',dmatcontents);
dmatcontents = cellfun(@(x)x(1:end-2),dmatcontents(ism),'uni',0);
funs = reshape([dmatcontents;repmat({[]},size(dmatcontents))],1,[]);
diagnostics = struct(funs{:});

rand('seed',0);
randn('seed',0);

p=1;

%% dmatref
diagnostics.dmatref = dmatref;
fprintf('%4i%%',p),p=p+2;

%% standardparset
diagnostics.standardparset = standardparset(0,.1);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% simuldiff
diagnostics.simuldiff = simuldiff(diagnostics.standardparset(1,:),250);
fprintf('\b\b\b\b\b%4i%%',p),p=p+3;

%% simuldiffwo
diagnostics.simuldiffwo = simuldiffwo(diagnostics.standardparset(2,:),...
    250,.025,[0 .2],.025,[.5 3]);
diagnostics.data = [[ones(250,1);ones(250,1)*2],...
    [diagnostics.simuldiff;diagnostics.simuldiffwo]];
fprintf('\b\b\b\b\b%4i%%',p),p=p+4;

%% multisimul
diagnostics.multisimul = multisimul;
fprintf('\b\b\b\b\b%4i%%',p),p=p+4;

%% multisimulwo
diagnostics.multisimulwo = multisimulwo(diagnostics.standardparset,20,...
    1,.05,[.05 .5],0.1,[.5 7]);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% isvaliddataset
[a b] = isvaliddataset(diagnostics.data);
diagnostics.isvaliddataset = {a b};
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% descriptives
diagnostics.descriptives = descriptives(diagnostics.data);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% splitdata
diagnostics.splitdata = splitdata(diagnostics.data);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% processdata
diagnostics.processdata = processdata(diagnostics.splitdata);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% multiestv4 (options)
diagnostics.multiestv4.options = multiestv4;
ns = 2;
o = orthpoly(ns,0);
e = orthpoly(ns,ns);
fv = nan(ns,9);
fv(ns,8:9) = [1 0];
diagnostics.multiestv4.options.DesignMatrix = {o o o e o o e e e};
diagnostics.multiestv4.options.nsets = ns;
diagnostics.multiestv4.options.npar = 9;
diagnostics.multiestv4.options.FixedValues = fv;
diagnostics.multiestv4.options.SpecificBias = repmat(.5,ns,1);
diagnostics.multiestv4.options.OutlierTreatment = 'b';
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% inpcheck
diagnostics.inpcheck = inpcheck(diagnostics.multiestv4.options);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% ewmav2
[co dmatewmaplot]= ewmav2(diagnostics.data,...
    diagnostics.inpcheck.EWMA.L,...
    diagnostics.inpcheck.EWMA.l,...
    diagnostics.inpcheck.EWMA.s);
diagnostics.ewmav2 = {co,dmatewmaplot};
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% outliertreatment
[data,oltrreport] = outliertreatment(diagnostics.data,...
    diagnostics.inpcheck);
diagnostics.outliertreatment = {data,oltrreport};
fprintf('\b\b\b\b\b%4i%%',p),p=p+3;

%% ezdiff
[v,a,Ter] = ezdiff(diagnostics.outliertreatment{1});
diagnostics.ezdiff = {v,a,Ter};
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% generateguess
diagnostics.generateguess = generateguess(diagnostics.inpcheck,...
    diagnostics.splitdata);
diagnostics.inpcheck.controls.large = diagnostics.generateguess.Guess;
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% smaller
diagnostics.smaller = smaller(diagnostics.inpcheck.controls);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% bigger
diagnostics.bigger = bigger(diagnostics.smaller);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% isconsistent
diagnostics.isconsistent = isconsistent(diagnostics.bigger);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% cdfdif
diagnostics.cdfdif = cdfdif(-.5:.5:3,0,diagnostics.standardparset(1,:));
fprintf('\b\b\b\b\b%4i%%',p),p=p+3;

%% invcdfdif
fprintf('  (please wait)')
diagnostics.invcdfdif = invcdfdif(diagnostics.cdfdif,0,...
    diagnostics.standardparset(1,:));
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%4i%%',p),p=p+12;

%% fitdiffv13
diagnostics.fitdiffv13 = fitdiffv13(diagnostics.standardparset(1,:),...
    diagnostics.processdata{1},1);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% multifitv4
diagnostics.multifitv4 = multifitv4(diagnostics.smaller.small,...
    diagnostics.processdata,diagnostics.smaller,1);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% iseye
diagnostics.iseye = {iseye(1),iseye(eye(3)),...
    iseye(eye(2,3)),iseye(0),iseye(NaN),iseye(ones(4))};
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% isgood
[a b] = isgood(diagnostics.standardparset);
[c d] = isgood(diagnostics.standardparset(:,1:6));
[e f] = isgood([.08 .01 .5 .075 .1 .01 -.2 .5 -.001]);
diagnostics.isgood = {a,b,c,d,e,f};
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% isnested
a = eye(3);
b = ones(3,1);
c = 0;
d = [0;1;1];
e = [b d];
f = [NaN;1;1];

in={isnested(a,a),isnested(a,b),isnested(a,c),...
    isnested(a,d),isnested(a,e),isnested(a,f)
    isnested(b,a),isnested(b,b),isnested(b,c),...
    isnested(b,d),isnested(b,e),isnested(b,f)
    isnested(c,a),isnested(c,b),isnested(c,c),...
    isnested(c,d),isnested(c,e),isnested(c,f)
    isnested(d,a),isnested(d,b),isnested(d,c),...
    isnested(d,d),isnested(d,e),isnested(d,f)
    isnested(e,a),isnested(e,b),isnested(e,c),...
    isnested(e,d),isnested(e,e),isnested(e,f)
    isnested(f,a),isnested(f,b),isnested(f,c),...
    isnested(f,d),isnested(f,e),isnested(f,f)};
% note the false positive in (2,4) by isnested(b,d)
diagnostics.isnested = in;
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% isposdef
diagnostics.isposdef = {isposdef(eye(4)),isposdef(-eye(4)),...
    isposdef([NaN 0;0 1]),isposdef(ones(2,1)),isposdef(zeros(1,2))};
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% chi2test
diagnostics.chi2test = chi2test([-1 0 NaN 10 10],[10 10 10 10 NaN]);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% multiestv4 (estimate)
a = diagnostics.inpcheck([1 1]);
a(1).OutlierTreatment = 'n';
a(1).Name = [mfilename '_1'];
a(2).Name = [mfilename '_2'];
fprintf('  (patience...)')
evalc('a = multiestv4(diagnostics.data,a);');
diagnostics.multiestv4.estimate = a;
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%4i%%',p),p=p+17;

%% fitlast
diagnostics.fitlast = fitlast(diagnostics.standardparset,...
    diagnostics.multiestv4.estimate);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% modelfittable
diagnostics.modelfittable = modelfittable(diagnostics.multiestv4.estimate);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% namepars
diagnostics.namepars = namepars(diagnostics.multiestv4.estimate);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% qtable
diagnostics.qtable = qtable(diagnostics.multiestv4.estimate);
fprintf('\b\b\b\b\b%4i%%',p),p=p+2;

%% Remove undesirables
diagnostics.multiestv4.estimate = [];
diagnostics.modelfittable = [];
diagnostics.cdfdif = [];
diagnostics.fitdiffv13 = [];
diagnostics.fitlast = [];
diagnostics.invcdfdif = [];
diagnostics.modelfittable = [];
diagnostics.multifitv4 = [];
diagnostics.namepars = [];
diagnostics.qtable = [];    
fprintf('\b\b\b\b\b%4i%%',p)

%% There! All done.
% In case we want to make a new benchmark file, this returns the results
% structure.
if nargin && strcmp(ret,'return')
    varargout{1} = diagnostics;
    fprintf('\b\b\b\b\bBenchmark data compiled.\n')
    return
end

%% Compare
c = load([dmatdir filesep 'test' filesep 'dmatdiagn.bmf'],'-mat');
v = isequalwithequalnans(diagnostics,c.dia);
d = [];
warning(stat)
fprintf('\b\b\b\b\b')

if ~v
    flnm = fieldnames(diagnostics);
    n = length(flnm);
    funs = reshape([flnm';repmat({false},1,n)],1,[]);
    d = struct(funs{:});
    r = zeros(1,n);
    for ctr=1:n
        d.(flnm{ctr}) = isequalwithequalnans(diagnostics.(flnm{ctr}),...
            c.dia.(flnm{ctr}));
        r(ctr) = d.(flnm{ctr});
    end
end

%% Report
fprintf(repmat('\b',1,59))
if ~nargout
    if v
        fprintf('The time is %s and all is well.',datestr(now,'HH:MM'))
    else
        fprintf(['DMAT performs without errors, but %i%% of its ',...
            'functions do not perform exactly as expected.\nIf you are ',...
            'experiencing problems, you can reinstall DMAT from the <a',...
            ' href="%s">DMAT website</a>.\n'],round(100*mean(~r)),dmsite)
        fprintf('Discrepancies occurred in:\n')
        disp(flnm(~r))
    end
else
    if nargout>0
        varargout{1} = v;
        if nargout>1
            varargout{2} = diagnostics;
            if nargout>2
                varargout{3} = d;
            end
        end
    end
end

fprintf('\nTEST_MAIN succesfully terminated after approximately %i seconds.\n',round(toc))