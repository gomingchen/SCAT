function En = energy(x,wintype,winamp,winlen)
%ENERGY   Short-time energy computation.
%   y = ENERGY(X,WINTYPE,WINAMP,WINLEN) computes the short-time enery of
%   the sequence X. 
%
%   WINTYPE defines the window type. RECTWIN, HAMMING, HANNING, and
%   BLACKAMN are the possible choices. WINAMP sets the amplitude of the
%   window and the length of the window is WINLEN. 
%   
%   See also RECTWIN, HAMMING, HANNING, BARTLETT, BLACKMAN.
%
%   Author: Nabin Sharma
%   Date: 2009/03/15

narginchk(1,4);

% generate the window
win = (winamp*(window(str2func(wintype),winlen))).';

% enery calculation
x2 = x.^2;
En = winconv(x2,wintype,win,winlen);
En = En(winlen/2:length(x)+winlen/2);       % truncate to relevant signal

% ensure only positive energy values
En(En<0) = 0;
