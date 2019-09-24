function [B,A] = fb_cheby1(N,R,Wn)
% FB_CHEBY1  creates a Chebychev type-I filterbank using specified input parameters
%
% The MxN output matrix contains M row vectors of time series data N samples long

for i = 1:size(Wn,1)
    [B(i,:),A(i,:)] = cheby1(N, R, Wn(i,:));
end
