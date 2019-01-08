function [ pOpt, R2, Range, hModel ] = FitNakaRushton( x, y , noise,SupplyGradient)
    % Written by Spero Nicholas and adapted by Peter J. Kohler, 2017
    % 
    % Description:	Fit Naka-Rushton function to Sweep data: 
    %               y = rMax / ( 1 + ( c50 / x )^exponent ) + baseline
    %               
    % 
    % Syntax:	[ pOpt, R2, Range, hModel ] = FitNakaRushton( x, y, SupplyGradient )
    % In:
    %   x               - 1 x n vector, containing the bin values 
    %    
    %   y               - n x m x q matrix, containing the sweep data, 
    %                    first dimension corresponds to the sweep bins, 
    %                    all other dimensions will be considered additional 
    %                    data sets and fit separately
    %   
    %   SupplyGradient - [true]/false, indicates whether or not to a
    %                    gradient should be supplied to fmincon
    % Out:
    % 	pOpt           - 4 x m x q x ... matrix, containing fit parameters, 
    %                   first dimension is [ c50, exponent, rMax, baseline ]
    %                   any additional dimensions correspond to additional data sets
    % 
    %   Range          - 1 x m x q matrix, the data range (min-max) for the fitted function,
    %                    same logic as pOpt.
    %
    %   R2             - 1 x m x q matrix, containing R2 goodness-of-fit estimates,
    %                    same logic as pOpt.
    %
    %   hModel         - handle of the model function used, will be the
    %                    same for all data sets
    
    if nargin < 3
        noise = zeros(size(y));
    else
    end
    
    if nargin < 4
        SupplyGradient = true;
    else
    end

    nx = size( x );
    if ~isvector( x ) || prod( nx ) == 1
        error( 'input x needs to be a vector' );
    end
    if any( x <= 0 )
        error( 'x values must be > 0' )
    end
    if nx(1) == 1
        x = x(:);
    end
    nx = numel( x );

    ny = size( y );
    if ny(1) ~= nx
        error( 'x & y length don''t agree' )
    end

    nPar = 4;
    pInit = zeros(nPar,1);					% [ c50; exponent; rMax; baseline ]
    pLB   = zeros(nPar,1);
    pUB   =   inf(nPar,1);
    pLB(3) = -Inf;						    % allow negative rMax
    if pLB(3) < 0
        Aie   = [ 0, 0, -1, -1 ];			% -rMax-baseline <= 0.   don't allow negative rMax+baseline
        bie   = 0;
    else
        Aie   = [];
        bie   = [];
    end
    xMin = min( x );
    xMax = max( x );

    pOpt = zeros( [ nPar , ny(2:end) ] );
    algs = { 'active-set', 'interior-point', 'sqp', 'trust-region-reflective' };
    opts = optimoptions( @fmincon, 'Display', 'notify', 'MaxFunEvals', 5e4 , 'MaxIter',500);
    switch SupplyGradient
        case 1		% supply gradient
            opts = optimoptions( opts, 'Algorithm', algs{3}, 'GradObj', 'on' );
        case 2		% don't
            opts = optimoptions( opts, 'Algorithm', algs{3}, 'GradObj', 'off' );
    end

    R2 = zeros( [ 1, ny(2:end) ] );
    Range = zeros( [ 1, ny(2:end) ] );

    nCurve = prod( ny(2:end) );
    % open parallel pool, if not already open
    if isempty(gcp)
        parpool;
    end
    % loop over the different datasets
    for iCurve = 1:nCurve
        if all(noise(:,iCurve) > 0)
            pLB(4) = min(noise(:,iCurve));
            pUB(4) = max(noise(:,iCurve));
        else
        end
        [ pOpt(:,iCurve), R2(iCurve), Range(iCurve) ] = doFit(x,y(:,iCurve),xMin,xMax,pLB,pUB,Aie,bie,nx,opts);
    end
    hModel = @(x,p) NRmodel(x,p);
end

function [pOpt,R2,Range] = doFit(x,y,xMin,xMax,pLB,pUB,Aie,bie,nx,opts)
        % set initial values
        SST = sum( ( y - mean( y, 1 ) ).^2 );
        % try full range of x
        [ SSEset(1), pInitSet(:,1) ] = setParsForXRange( x, y, SST, xMin, xMax, pLB );
        % right half of x range
        [ SSEset(2), pInitSet(:,2) ] = setParsForXRange( x, y, SST, (xMin + xMax ) / 2, xMax, pLB );
        % left half of x range
        [ SSEset(3), pInitSet(:,3) ] = setParsForXRange( x, y, SST, xMin, ( xMin + xMax ) / 2, pLB );
        % middle half of xrange
        [ SSEset(4), pInitSet(:,4) ] = setParsForXRange( x, y, SST, 0.75 * xMin + 0.25 * xMax, 0.25 * xMin + 0.75 * xMax, pLB );
        % find best value, as smallest sum-of-squared errors
        [~, bestIdx ] = min(SSEset);
        pInit = pInitSet(:,bestIdx);
        % optimize
        pOpt = fmincon( @(p)NRobjFcn(p,x,y,opts), pInit, Aie, bie, [], [], pLB, pUB, @(p)NRconFcn(p,xMin,xMax,nx), opts );
        if pOpt(3) == 0
            % if rMax == 0, c50 & exponent are meaningless
            pOpt(1) = 0;		% set c50 to middle of range?
            pOpt(2) = 0;
        end
        % get actual excursion of data range
        Range = NRmodel( xMax, pOpt ) - NRmodel( xMin, pOpt );
        % goodness of fit: 1-SSE/SST
        R2 = 1 - sum( ( y - NRmodel( x, pOpt ) ).^2 ) / SST;
end

function [ SSEinit,pInit ] = setParsForXRange( x, y, SST, xLow, xHigh, pLB )
    Fcrit = [ 0.01; 0.99 ];			% actual criterion range of delta-y
    Fcrit(:) = 1./Fcrit - 1;		% transform to values used in equation
    % set c50 & exponent so NRmodel( [xLow;xHigh], [c50;exp;1;0] ) == Fcrit
    pInit(2) = log( Fcrit(1) / Fcrit(2) ) / log( xHigh / xLow );
    pInit(1) = xLow * Fcrit(1)^( 1 / pInit(2) );
    % set rMax & baseline to least squares data fit given c50 & exponent
    if pLB(3) < 0
        pInit(3:4) = [ NRmodel( x, [ pInit(1:2)'; 1; 0 ] ), ones(size(x)) ] \ y;
    else
        pInit(3:4) = lsqnonneg( [ NRmodel( x, [ pInit(1:2)'; 1; 0 ] ), ones(size(x)) ], y, optimset( 'Display', 'notify' ) );
    end
    SSEinit = sum( ( y - NRmodel( x, pInit ) ).^2 );
end

function y = NRmodel( x, p )
    y = p(3) ./ ( ( p(1) ./ x ).^p(2) + 1 ) + p(4);		% x > 0
end

function [ f, g ] = NRobjFcn( p, x, y, opts )
    r = y - NRmodel( x, p );
    f = sum( r.^2 );
    if strcmp( opts.GradObj, 'on' )
        g = [ 
            sum(  2 * r * p(3) ./ ( 1 + (p(1)./x).^p(2) ).^2 * p(2) .* (p(1)./x).^( p(2) - 1 ) ./ x )
            sum(  2 * r * p(3) ./ ( 1 + (p(1)./x).^p(2) ).^2 .* (p(1)./x).^p(2) .* log(p(1)./x) )
            sum( -2 * r        ./ ( 1 + (p(1)./x).^p(2) ) )
            sum( -2 * r )
        ];
    end
end

function [ cie, ceq ] = NRconFcn( p, xMin, xMax, nx )
    % cie(p) <= 0, ceq(p) == 0	
    cie = [ 
        NRmodel( xMin, [ p(1:2); 1; 0 ] ) - 0.10															%  lowest bin invisible
        0.10 - NRmodel( xMax, [ p(1:2); 1; 0 ] )															% highest bin not invisible.  zero it w/ p(3) if needed
        ( xMax - xMin ) / ( nx - 1 ) - ( (1/(1/0.99-1))^(1/p(2)) - (1/(1/0.01-1))^(1/p(2)) ) * p(1)			% not too steep
        ];
    ceq = [];
end