function [forward, feedback] = fb_gammatone(fs, fc, ERB)
% FB_GAMMATONE  creates a gammatone filterbank using specified input parameters
%
% The MxN output matrix contains M row vectors of time series data N samples long
%
% Note:  Following equations are taken from Malcolm Slaney's GammaToneMake.m, 
%   which is from Apple TR #35:
% "An Efficient Implementation of the Patterson-Holdsworth Cochlear Filter Bank."

T = 1/fs;
B=1.019*2*pi*ERB; %in rad here - note to self: some models require B in Hz (NC)
gain = abs((-2*exp(4*1i*fc*pi*T)*T + 2*exp(-(B*T) + 2*1i*fc*pi*T).*T.*(cos(2*fc*pi*T) - sqrt(3 - 2^(3/2))*sin(2*fc*pi*T))) .*(-2*exp(4*i*fc*pi*T)*T +2*exp(-(B*T) + 2*i*fc*pi*T).*T.*(cos(2*fc*pi*T) + sqrt(3 - 2^(3/2)) *sin(2*fc*pi*T))).*(-2*exp(4*i*fc*pi*T)*T +2*exp(-(B*T) + 2*i*fc*pi*T).*T.*(cos(2*fc*pi*T) -sqrt(3 + 2^(3/2))*sin(2*fc*pi*T))) .*(-2*exp(4*i*fc*pi*T)*T+2*exp(-(B*T) + 2*i*fc*pi*T).*T.*(cos(2*fc*pi*T) + sqrt(3 + 2^(3/2))*sin(2*fc*pi*T))) ./(-2 ./ exp(2*B*T) - 2*exp(4*i*fc*pi*T) +2*(1 + exp(4*i*fc*pi*T))./exp(B*T)).^4);
feedback=zeros(length(fc),9);
forward=zeros(length(fc),5);
forward(:,1) = T^4 ./ gain;
forward(:,2) = -4*T^4*cos(2*fc*pi*T)./exp(B*T)./gain;
forward(:,3) = 6*T^4*cos(4*fc*pi*T)./exp(2*B*T)./gain;
forward(:,4) = -4*T^4*cos(6*fc*pi*T)./exp(3*B*T)./gain;
forward(:,5) = T^4*cos(8*fc*pi*T)./exp(4*B*T)./gain;
feedback(:,1) = ones(length(fc),1);
feedback(:,2) = -8*cos(2*fc*pi*T)./exp(B*T);
feedback(:,3) = 4*(4 + 3*cos(4*fc*pi*T))./exp(2*B*T);
feedback(:,4) = -8*(6*cos(2*fc*pi*T) + cos(6*fc*pi*T))./exp(3*B*T);
feedback(:,5) = 2*(18 + 16*cos(4*fc*pi*T) + cos(8*fc*pi*T))./exp(4*B*T);
feedback(:,6) = -8*(6*cos(2*fc*pi*T) + cos(6*fc*pi*T))./exp(5*B*T);
feedback(:,7) = 4*(4 + 3*cos(4*fc*pi*T))./exp(6*B*T);
feedback(:,8) = -8*cos(2*fc*pi*T)./exp(7*B*T);
feedback(:,9) = exp(-8*B*T);
