%% SOME_PARAMETERS  Define some variables for use in DMAT (demo script)

%% First, define a parameter set from the STANDARDPARSET function
param = standardparset(1);
disp(param)

%% Now, get its dimensions
[nsets npar] = size(param);

%% Make a matrix of fixed values
% If we don't want to fix any parameters to certain values, keep all
% elements of this matrix NaN. Suppose we want to fix the first and second
% conditions' drift rates to zero.
fixvals = nan(nsets,npar);
fixvals(1:2,7)=0;
disp(fixvals)

%% Make a vector of bias values
% We don't want any biases.
bias = nan(nsets,1);

%% Now, create a design matrix
% The design matrix is contained in a cell matrix, wherein each cell is the
% design matrix for a certain parameter.
designmat = {orthpoly(nsets,0);  % boundary separation
    orthpoly(nsets,0);           % mean nondecision time
    orthpoly(nsets,0);           % trial to trial stdev of drift rate
    orthpoly(nsets,0);           % mean starting point
    orthpoly(nsets,0);           % trial to trial range of starting point
    orthpoly(nsets,0);           % trial to trial range of nondecision time
    orthpoly(nsets,0:1)};        % drift rate
%%
% Notice the |orthpoly| function, which builds matrices of orthogonal
% polygons, for a given number of conditions (nsets) and a given number of
% degrees (here, the 0th and 1st degree). Also notice that this introduces
% a possible conflict between the fixed values matrix and the design matrix
% (both are imposing restrictions).

%% The CONTROLS structure contains all the most important information
% When building a structure, don't forget to put the curly braces around
% cells (look at how |designmat| is included). Otherwise a structure
% array will be created.
controls = struct('large',param,...
    'fixvals',fixvals,...
    'bias',bias,...
    'designmat',{designmat});
controls = smaller(controls);
disp(controls)
%%
% Notice here how the conflict between fixed values and the design was
% resolved by the |smaller| function (it overwrote the |fixvals| field).
disp(controls.fixvals)

%% Quick clean-up
% Just to see that the variables introduced to the workspace are consistent
% (this is not strictly necessary, just clean and to avoid confusion):
fixvals(:) = NaN;

%% And build an OPTIONS structure, ready for DMAT
options = struct(...
    'DesignMatrix',{designmat},...
    'FixedValues',fixvals,...
    'SpecificBias',bias,...
    'Name','Linear effect on drift rate');
disp(options)

%%
%%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.