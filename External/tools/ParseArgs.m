function varargout = ParseArgs(vargin,varargin)
% ParseArgs
% 
% Description:	parse a varargin cell of optional arguments and options.
%				optional arguments are arguments that do not need to be included
%				in a function call and assume a default value if omitted.
%				options are 'key'/value pairs that come at the end of a function
%				call argument list.
% 
% Syntax:	[v1,v2,...,vn,[opt]] = ParseArgs(vargin,d1,...,dn[,opt1,opt1def,...,optM,optMdef])
%
% In:
%	vargin		- the varargin cell
%	dK			- the default value of the Kth optional argument
%	[optJ]		- the name of the Jth option
%	[optJdef]	- the default value of the Jth option
% 
% Out:
%	vK			- the value of the Kth optional argument
%	[opt]		- a struct of option values. options specified by the user but
%				  not given default values are place in opt.opt_extra
%
% Note:	if the user calls the function like this:
%		func(v1,...,vN-1,vN,opt1,val1,...,optM,valM), vN-1 might possibly have
%		the same value as one of the option names, and that option wasn't
%		explicitly set in the options section, then vN-1 will be confused with
%		the option.
% 
% Updated: 2015-12-10
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
nDefault	= numel(varargin);
nOut		= nargout;
bOption		= nDefault~=nOut;

nUser		= numel(vargin);
nArgument	= nOut - bOption;

%reshape just to make sure
	vargin		= reshape(vargin,nUser,1);
	varargin	= reshape(varargin,nDefault,1);

if bOption
	varargout	= DoParseOpt;
else
	varargout	= DoParse;
end

%------------------------------------------------------------------------------%
function out = DoParseOpt
	%split the input between optional arguments and options
		%split the defaults between optional arguments and options
			defArgument	= varargin(1:nArgument);
			defOptKey	= varargin(nArgument+1:2:nDefault);
			defOptVal	= varargin(nArgument+2:2:nDefault);
			nDefOpt		= numel(defOptKey);
		%split user input between optional arguments and options
			if nUser==0
			%user didn't specify any arguments, just use the defaults
				if isempty(defOptVal)
					opt	= struct;
				else
					opt	= cell2struct(defOptVal,defOptKey);
				end
				
				opt.isoptstruct	= true;
				opt.opt_extra	= struct;
				
				out	= [defArgument; opt];
				
				return;
			else
			%the user specified arguments explicitly
				%head backward to find the earliest user argument that might be
				%an option key
					bOptionsSpecified	= false;
					for kChar=nUser-1:-2:1
						bChar	= ischar(vargin{kChar});
						if bChar
							bOptionsSpecified	= true;
						else
							kChar	= kChar+2;
							break;
						end
					end
				%what do we have?
					if bOptionsSpecified
						kKeyUserPossible	= kChar:2:nUser-1;
						userOptKeyPossible	= vargin(kKeyUserPossible);
						userOptValPossible	= vargin(kKeyUserPossible+1);
						nUserOptPossible	= numel(userOptKeyPossible);
						
						%which strings match default options?
							bOptKeyMatch	= false(nUserOptPossible,1);
							kDefOptKey		= zeros(nUserOptPossible,1);
							
							for kU=1:nUserOptPossible
								for kD=1:nDefOpt
									if strcmp(userOptKeyPossible{kU},defOptKey{kD})
										bOptKeyMatch(kU)	= true;
										kDefOptKey(kU)		= kD;
										break;
									end
								end
							end
						
						%the first user opt key is either:
						%	a) the first possible key that matches a default key
						%	b) the first possible key after the number of
						%	   optional arguments,
						%whichever comes earlier.
						%
						%in case of b, shift back by one key/value pair if that
						%would leave dangling values in between the last
						%optional argument and the first option key
						
						%position in user vargin of first matching option key
							kFirstOptKeyUserMatch	= kChar + 2*(find(bOptKeyMatch,1)-1);
							if isempty(kFirstOptKeyUserMatch)
								kFirstOptKeyUserMatch	= inf;
							end
						
						%position of first option key based on optional
						%arguments
							kFirstOptKeyFromArgument	= max(kChar,nArgument+1);
						
						if kFirstOptKeyUserMatch <= kFirstOptKeyFromArgument
						%case a above
							kFirstOptKey	= kFirstOptKeyUserMatch;
						else
						%case b above
							if kFirstOptKeyFromArgument>nUser
							%no options specified
								kFirstOptKey	= nUser+1;
							elseif ~ischar(vargin{kFirstOptKeyFromArgument})
							%we landed on an opt value rather than an opt key.
							%push one index backward.
								kFirstOptKey	= kFirstOptKeyFromArgument - 1;
							else
								kFirstOptKey	= kFirstOptKeyFromArgument;
							end
						end
						
						%separate the matching and extra options
							%first just get the actual option key/value pairs
								for kOptFirst=1:nUserOptPossible
									if kKeyUserPossible(kOptFirst)>=kFirstOptKey
										break;
									end
								end
								
								%kOptFirst	= find(kKeyUserPossible>=kFirstOptKey,1);
								kOptKeep	= kOptFirst:nUserOptPossible;
								nUserOpt	= nUserOptPossible - kOptFirst + 1;
								
								userOptKey		= userOptKeyPossible(kOptKeep);
								userOptVal		= userOptValPossible(kOptKeep);
								bOptKeyMatch	= bOptKeyMatch(kOptKeep);
								kDefOptKey		= kDefOptKey(kOptKeep);
							
							%now extract the ones that don't match
								bOptKeyNoMatch	= ~bOptKeyMatch;
								
								bExtraOptions	= any(bOptKeyNoMatch);
								
								if bExtraOptions
									extraOptKey	= userOptKey(bOptKeyNoMatch);
									extraOptVal	= userOptVal(bOptKeyNoMatch);
									
									userOptKey(bOptKeyNoMatch)	= [];
									userOptVal(bOptKeyNoMatch)	= [];
									kDefOptKey(bOptKeyNoMatch)	= [];
									
									nUserOpt	= numel(userOptKey);
								end
						
						%eliminate explicitly unspecified options
							bUserUnspecified	= false(nUserOpt,1);
							for kU=1:nUserOpt
								bUserUnspecified(kU)	= isempty(userOptVal{kU});
							end
							
							userOptKey(bUserUnspecified)	= [];
							userOptVal(bUserUnspecified)	= [];
							kDefOptKey(bUserUnspecified)	= [];
							nUserOpt						= numel(userOptKey);
						
						%replace default options with user specified options
							optKey	= defOptKey;
							optVal	= defOptVal;
							
							for kO=1:nUserOpt
								optVal{kDefOptKey(kO)}	= userOptVal{kO};
							end
					else
					%no options were specified
						kFirstOptKey	= nUser+1;
						
						optKey	= defOptKey;
						optVal	= defOptVal;
						
						bExtraOptions	= false;
					end
				
				%user optional arguments
					nUserArgument	= min(nArgument,kFirstOptKey-1);
					userArgument	= vargin(1:nUserArgument);
			end
	
	%parse optional arguments
		out	= defArgument;
		
		for kU=1:nUserArgument
			if ~isempty(userArgument{kU})
				out{kU}	= userArgument{kU};
			end
		end
	%parse the options
		%get the extra options specified by the user
			if bExtraOptions
				try
				%first assume there are no duplicate field names
					opt_extra	= cell2struct(extraOptVal, extraOptKey);
				catch me
					switch me.identifier
						case 'MATLAB:DuplicateFieldName'
						%ok, this will be a bit slower
							[extraOptKey,kUnique]	= unique(extraOptKey);
							extraOptVal				= extraOptVal(kUnique);
							opt_extra				= cell2struct(extraOptVal, extraOptKey);
						otherwise
							disp(errID);
							rethrow(me);
					end
				end
			else
				opt_extra	= struct;
			end
		
		%construct the options struct
			if ~isempty(optVal)
				opt	= cell2struct(optVal,optKey);
			else
				opt	= struct;
			end
		
		%mark opt as an options struct (see optstruct)
			opt.isoptstruct	= true;
			opt.opt_extra	= opt_extra;
		
		out{end+1}	= opt;
end
%------------------------------------------------------------------------------%
function out = DoParse
	out	= vargin;
	if nUser<nArgument
		out{nArgument}	= [];
	end
	
	for k=1:nArgument
		if isempty(out{k})
			out{k}	= varargin{k};
		end
	end
end
%------------------------------------------------------------------------------%

end
