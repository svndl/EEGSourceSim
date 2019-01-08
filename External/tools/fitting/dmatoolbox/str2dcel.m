function [dcel dcel2] = str2dcel(str)
%STR2DCEL  Convert data structure into DCEL format
%  [DCEL DCEL2] = STR2DCEL(STRUCT), where STRUCT has all the fields
%            .CorrectObs          .Number of observations in each X=1 bin
%            .IncorrectObs        .Number of observations in each X=0 bin
%            .CorrectEdges        .Edges of X=1 time bins (in sec)
%            .IncorrectEdges      .Edges of X=0 time bins (in sec)
%            .Lower               .Minimal RT (optional)
%            .Upper               .Maximal RT (optional)
%  returns DCEL and DCEL2.
%
%  See also DCEL2STR, QUANTEST.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

n = numel(str);
dcel2 = cell(n,1);
dcel = dcel2;

for ctr = 1:n
    if isfield(str(ctr),'Lower')&&isfield(str(ctr),'Upper')
        lo = str(ctr).Lower;
        up = str(ctr).Upper;
        limits = [lo up];
    else
        lo = .95*min([str(ctr).CorrectEdges(:);str(ctr).IncorrectEdges(:)]);
        up = 3*max([str(ctr).CorrectEdges(:);str(ctr).IncorrectEdges(:)]);
        limits = [NaN NaN];
    end
    sum1 = sum(str(ctr).CorrectObs);
    sum0 = sum(str(ctr).IncorrectObs);
    dcel2{ctr} = {str(ctr).CorrectEdges(:)'
        str(ctr).IncorrectEdges(:)'
        str(ctr).CorrectObs(:)'
        str(ctr).IncorrectObs(:)'
        sum1(:)'
        sum0(:)'
        sum1+sum0
        limits(:)'};
    nbins1 = length(str(ctr).CorrectEdges(:));
    bins1 = [lo(:); str(ctr).CorrectEdges(:); up(:)];
    nbins0 = length(str(ctr).IncorrectEdges(:));
    bins0 = [lo(:); str(ctr).IncorrectEdges(:); up(:)];
    obs1 = [0; cumsum(str(ctr).CorrectObs(:))];
    obs0 = [0; cumsum(str(ctr).IncorrectObs(:))];
    cors = [];
    errs = [];
    for a=1:nbins1+1
        sub = bins1(a)*rand(1,obs1(a+1)-obs1(a))+lo;
        cors = [cors;sub(:)];
    end
    for a=1:nbins0+1
        sub = bins0(a)*rand(1,obs0(a+1)-obs0(a))+lo;
        errs = [errs;sub(:)];
    end
    dcel{ctr} = [ones(length(cors),1) cors(:)
        zeros(length(errs),1) errs(:)];
end