function s = pairs2struct(varargin)
%PAIRS2STRUCT  Turns field-value pairs into a structure
%   S = PAIRS2STRUCT('field1','value1','field2','value2',...), returns
%   S.(field1)=value1, etc.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

if rem(nargin,2)
   error('DMAT:pairs2struct:badInput','Inputs need to come in pairs.')
end

names = fieldnames(multiestv4);
s = struct;

for ctr = 1:2:nargin-1
   newname = varargin{ctr};
   if ~ismember(newname,names)
      error('DMAT:pairs2struct:badInput',...
         'Unknown input field ''%s''.',newname)
   end
   s.(newname) = varargin{ctr+1};
end
