function str = dcel2str(dcel2)
%DCEL2STR  Convert DCEL format into data structure
%  STRUCT = DCEL2STR(DCEL2), returns STRUCT with all the fields
%            .CorrectObs          .Number of observations in each X=1 bin
%            .IncorrectObs        .Number of observations in each X=0 bin
%            .CorrectEdges        .Edges of X=1 time bins (in sec)
%            .IncorrectEdges      .Edges of X=0 time bins (in sec)
%            .Lower               .Minimal RT (optional)
%            .Upper               .Maximal RT (optional)
%
%  See also STR2DCEL, QUANTEST.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

n = length(dcel2);

for ctr = 1:n
    str(ctr).CorrectEdges = dcel2{ctr}{1};
    str(ctr).IncorrectEdges = dcel2{ctr}{2};
    str(ctr).CorrectObs = dcel2{ctr}{3};
    str(ctr).IncorrectObs = dcel2{ctr}{4};
    str(ctr).Lower = dcel2{ctr}{8}(1);
    str(ctr).Upper = dcel2{ctr}{8}(2);
end