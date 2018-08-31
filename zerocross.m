function [f,ndx] = zerocross(vector)
%ZEROCROSS Finds number of zerocrosses
%    Zerocross simply reports the number of times
%    the input vector crosses the zero boundary. If two
%    outputs defined, also returns indices
%    (the latter index of crossing).
%
%    Usage
%      [NUMBER,INDEX]=ZEROCROSS(X);

% (c) mikko nissi 2012, <nissi@cmrr.umn.edu>

ndx=find(abs(diff(sign(vector))));
if~isempty(ndx)
    ndx=ndx+1;
    f=length(ndx);
else
    return
end


