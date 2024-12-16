function depth = treedpthX(t)
%TREEDPTH Tree depth. 
%   D = treedpthX(T) returns the depth D of the tree T.
%
%   See also WTREEMGR.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 14-May-2003.
%   Copyright 1995-2004 The MathWorks, Inc.
% $Revision: 1.1 $

depth = wtreemgr('depth',t);