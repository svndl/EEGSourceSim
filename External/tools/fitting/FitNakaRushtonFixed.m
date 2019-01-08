function [ pOpt, R2, hModel ] = FitNakaRushtonEq( x, y, yBaseline )
% Fit multi-dimensional data array "y" vs. vector "x" with Naka Rushton equation
%
% F( x | c50, exponent, rMax, baseline ) = rMax ./ ( 1 + ( c50 ./ x).^exponent ) + baseline
%
% USAGE:
% [ pOpt, R2, hModel ] = FitNakaRushtonEq( x, y, baseline )
%
% Inputs:
% x = vector
% y = array of any dimensions as long as size(y,1) = numel(x)
%
% Outputs:
% pOpt   = optimal parameters.  array matching dimesions of y, except size(pOpt,1) = 5
%          columns of pOpt = [ c50; exponent; rMax; baseline; F(xMax)-F(xMin) ]
% R2     = coefficient of determination for each fit
% hModel = function handle for model.  yFit = hModel(x,p)

% note: parfor won't work the way this is written

% Validate Inputs
narginchk( 2, 3 )
if nargin == 2
	yBaseline = [];
end
% -- x
nx = size( x );
if ~isvector( x ) || prod( nx ) == 1
	error( 'input x needs to be a vector' );
end
if any( x < 0 )
	error( 'x values must be >= 0' )
end
if nx(1) == 1
	x = x(:);
end
nx = numel( x );
% -- y
ny = size( y );
if ny(1) ~= nx
	error( 'x & y length don''t agree' )
end
xMin = min( x );
xMax = max( x );

% Initialize parameters, set bounds & inequality contstraints & optimization controls
if isempty( yBaseline )
	nPar = 4;								% #parameters that are optimized
else
	nPar = 3;								% user-supplied baseline
end
p4    = 0;
pInit = zeros(nPar,1);					% [ c50; exponent; rMax; baseline ].  y = rMax / ( 1 + ( c50 / x )^exponent ) + baseline
pLB   = zeros(nPar,1);
pUB   =   inf(nPar,1);
pLB(3) = -Inf;								% allow negative rMax

	% keep c50 in data range
	pLB(1) = xMin;
	pUB(1) = xMax;

if pLB(3) < 0
	if nPar == 4
		Aie = [ 0, 0, -1, -1 ];			% -rMax-baseline <= 0.   don't allow negative rMax+baseline
		bie = 0;
	else
		Aie = [ 0, 0, -1 ];
		bie = p4;
	end
else
	Aie = [];
	bie = [];
end
pOpt = zeros( [ nPar + 1, ny(2:end) ] );
algs = { 'active-set', 'interior-point', 'sqp', 'trust-region-reflective' };
opts = optimoptions( @fmincon, 'Display', 'notify', 'MaxFunEvals', 5e3 );
switch 1
case 1		% supply gradient
	opts = optimoptions( opts, 'Algorithm', algs{3}, 'GradObj', 'on' );				% can't use opts.UseParallel=true unless opts.SpecifyObjectiveGradient=false or opts.SpecifyContstraintGradient=false
case 2		% don't
	opts = optimoptions( opts, 'Algorithm', algs{3}, 'GradObj', 'off' );
end
R2 = zeros( [ 1, ny(2:end) ] );

% if strcmp( opts.GradObj, 'on' ) && any( x == 0 )
% 	warning( 'can''t supply analytic gradient when x = 0' )
% 	opts = optimoptions( opts, 'GradObj', 'off' );
% end

% Run optimizations
nCurve = prod( ny(2:end) );
for iCurve = 1:nCurve

	if nPar == 3
		p4(:) = yBaseline(iCurve);
		if ~isempty( bie )
			bie(:) = p4;
		end
	end

	% set initial values
	SST = sum( ( y(:,iCurve) - mean( y(:,iCurve), 1 ) ).^2 );
	setInitPars
	% optimize
	pOpt(1:nPar,iCurve) = fmincon( @NRobjFcn, pInit, Aie, bie, [], [], pLB, pUB, @NRconFcn, opts );
	if pOpt(3,iCurve) == 0
		% if rMax == 0, c50 & exponent are meaningless
		pOpt(1,iCurve) = 0;		% set c50 to middle of range?
		pOpt(2,iCurve) = 0;
	end
	% get model curve excursion over data range
	pOpt(nPar+1,iCurve) = NRmodel( xMax, pOpt(1:nPar,iCurve) ) - NRmodel( xMin, pOpt(1:nPar,iCurve) );
	% goodness of fit - coefficient of determination
	R2(iCurve) = 1 - sum( ( y(:,iCurve) - NRmodel( x, pOpt(1:nPar,iCurve) ) ).^2 ) / SST;
	
end

hModel = @(x,p) NRmodel(x,p);

return


	function setInitPars
		Fcrit = [ 0.01; 0.99 ];			% actual criterion range of delta-y
		Fcrit(:) = 1./Fcrit - 1;		% transform to values used in equation
		fInit = Inf;
		% try full range of x
		setParsForXRange( xMin, xMax );
		fBest = fInit;
		pBest  = pInit;
		% right half of x range
		setParsForXRange( ( xMin + xMax ) / 2, xMax );
		if fInit < fBest
			fBest(:) = fInit;
			pBest(:)  = pInit;
		end
		% left half of x range
		setParsForXRange( xMin, ( xMin + xMax ) / 2 );
		if fInit < fBest
			fBest(:) = fInit;
			pBest(:)  = pInit;
		end
		% middle half of xrange
		setParsForXRange( 0.75 * xMin + 0.25 * xMax, 0.25 * xMin + 0.75 * xMax );
		if fInit > fBest
			pInit(:) = pBest;
		end
		return
		function setParsForXRange( xLow, xHigh )
			if xLow == 0
				pInit(1) = xHigh / 2;
				pInit(2) = log( Fcrit(2) ) / log( 0.5 );
			else
				% set c50 & exponent so NRmodel( [xLow;xHigh], [c50;exp;1;0] ) == Fcrit
				pInit(2) = log( Fcrit(1) / Fcrit(2) ) / log( xHigh / xLow );
				pInit(1) = xLow * Fcrit(1)^( 1 / pInit(2) );
			end
			% set rMax & baseline to least squares data fit given c50 & exponent
			if nPar == 4
				if pLB(3) < 0
					pInit(3:4) = [ NRmodel( x, [ pInit(1:2); 1; 0 ] ), ones(nx,1) ] \ y(:,iCurve);
				else
					pInit(3:4) = lsqnonneg( [ NRmodel( x, [ pInit(1:2); 1; 0 ] ), ones(nx,1) ], y(:,iCurve), optimset( 'Display', 'notify' ) );
				end
			else
				if pLB(3) < 0
					pInit(3) = NRmodel( x, [ pInit(1:2); 1 ] ) \ y(:,iCurve);
				else
					pInit(3) = lsqnonneg( NRmodel( x, [ pInit(1:2); 1 ] ), y(:,iCurve), optimset( 'Display', 'notify' ) );
				end
			end
			fInit(:) = sum( ( y(:,iCurve) - NRmodel( x, pInit ) ).^2 );
		end
	end
	
	function y = NRmodel( x, p )
		if numel( p ) == 4
			y0 = p(4);
		else
			y0 = p4;
		end
		y = p(3) ./ ( ( p(1) ./ x ).^p(2) + 1 ) + y0;		% x > 0
		kx0 = x == 0;
		if any( kx0 )
			if p(2) > 0
				y(kx0) = y0;
			elseif p(2) < 0						% negative exponent shouldn't ever happen unless fmincon tries out-of-bounds parameters
				y(kx0) = p(3) + y0;
			else
				y(kx0) = p(3)/2 + y0;
			end
		end
	end

	function [ f, g ] = NRobjFcn( p )
		r = y(:,iCurve) - NRmodel( x, p );
		f = sum( r.^2 );
		if strcmp( opts.GradObj, 'on' )
			kx0 = x == 0;
			if any( kx0 )
				if p(2) > 0
					g = [ 
						sum(  2 * r(~kx0) * p(3) ./ ( 1 + (p(1)./x(~kx0)).^p(2) ).^2 * p(2) .* (p(1)./x(~kx0)).^( p(2) - 1 ) ./ x(~kx0) )
						sum(  2 * r(~kx0) * p(3) ./ ( 1 + (p(1)./x(~kx0)).^p(2) ).^2 .* (p(1)./x(~kx0)).^p(2) .* log(p(1)./x(~kx0)) )
						sum( -2 * r(~kx0)        ./ ( 1 + (p(1)./x(~kx0)).^p(2) ) )
						sum( -2 * r(~kx0) ) - 2*r(kx0)
					];
				elseif p(2) < 0						% negative exponent shouldn't ever happen unless fmincon tries out-of-bounds parameters
					g = [ 
						sum(  2 * r(~kx0) * p(3) ./ ( 1 + (p(1)./x(~kx0)).^p(2) ).^2 * p(2) .* (p(1)./x(~kx0)).^( p(2) - 1 ) ./ x(~kx0) )
						sum(  2 * r(~kx0) * p(3) ./ ( 1 + (p(1)./x(~kx0)).^p(2) ).^2 .* (p(1)./x(~kx0)).^p(2) .* log(p(1)./x(~kx0)) )
						sum( -2 * r(~kx0)        ./ ( 1 + (p(1)./x(~kx0)).^p(2) ) ) - 2*r(kx0)
						sum( -2 * r(~kx0) ) - 2*r(kx0)
					];
				else
					g = [ 
						sum(  2 * r(~kx0) * p(3) ./ ( 1 + (p(1)./x(~kx0)).^p(2) ).^2 * p(2) .* (p(1)./x(~kx0)).^( p(2) - 1 ) ./ x(~kx0) )
						sum(  2 * r(~kx0) * p(3) ./ ( 1 + (p(1)./x(~kx0)).^p(2) ).^2 .* (p(1)./x(~kx0)).^p(2) .* log(p(1)./x(~kx0)) )
						sum( -2 * r(~kx0)        ./ ( 1 + (p(1)./x(~kx0)).^p(2) ) ) - r(kx0)
						sum( -2 * r(~kx0) ) - 2*r(kx0)
					];
				end
			else
				g = [ 
					sum(  2 * r * p(3) ./ ( 1 + (p(1)./x).^p(2) ).^2 * p(2) .* (p(1)./x).^( p(2) - 1 ) ./ x )
					sum(  2 * r * p(3) ./ ( 1 + (p(1)./x).^p(2) ).^2 .* (p(1)./x).^p(2) .* log(p(1)./x) )
					sum( -2 * r        ./ ( 1 + (p(1)./x).^p(2) ) )
					sum( -2 * r )
				];
			end
			if nPar == 3
				g = g(1:3);
			end
		end
	end

	function [ cie, ceq ] = NRconFcn( p )
		% cie(p) <= 0, ceq(p) == 0
		fVisible = 0.1;
		fMaxRise = 0.99;
		fMinRise = 1 - fMaxRise;
		minBinRise = 2;		% fewest bins (assuming even spacing) allowed for rise from fMinRise to fMaxRise
		cie = [ 
			NRmodel( xMin, [ p(1:2); 1; 0 ] ) - fVisible																				%  lowest bin invisible
			fVisible - NRmodel( xMax, [ p(1:2); 1; 0 ] )																				% highest bin not invisible.  zero it w/ p(3) if needed
			minBinRise * ( xMax - xMin ) / ( nx - 1 ) - ( 1/(1/fMaxRise-1)^(1/p(2)) - 1/(1/fMinRise-1)^(1/p(2)) ) * p(1)			% not too steep
			];
		ceq = [];
	end

end