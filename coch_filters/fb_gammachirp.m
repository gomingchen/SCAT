function [B,A] = fb_gammachirp(fs, fc, ERB, fMax)
% FB_GAMMACHIRP  creates a gammachirp filterbank using specified input
% parameters
%
% Irino & Patterson use the Power Spectrum Model of Masking (Fletcher,
% 1940; Patterson, 1976) using experimental data to determine the gammatone
% filter shape.  Here, we attempt to use Prony's method (L. B. Jackson,
% 1996, Ch. 10) to estimate an ARMA filter given the impulse response,
% input signal amplitude, and default Gammatone parameters for a bat's
% cochlea model.
%
% The MxN output matrix contains M row vectors of time series data N
% samples long

%P_s = 40;       % SPL (dB)

% default gammachirp parameters
a = 1;           % impulse response amplitude
b = 1.68;        % distribution parameter
n = 4;           % distribution parameter
%c = 3.38 - 0.107*P_s   % FM parameter
c = 1;  % for gammatone FB
phi = 0;          % initial phase

%fc = 100000;       % asymptotic frequency
%ERB = 24.7*(1+0.00437*fc);      % equivalent rectangular bandwidth
%generate the complex gammachirp impulse response
%t = [1:400]./(20*fc(i));    % always keep 20 periods of fc

t = 1:round(fMax/20)/fs;
for i=1:length(fc)
    g_c = a*t.^(n-1) .* exp(-2*pi*b*ERB(i)*t) .* exp(1i*c*log(t) + 1i*phi);  % envelope
    g_c2 = exp(1i*2*pi*fc(i)*t) .* g_c; % frequency carrier
    
    % use prony to generate TF model
tic;    [B(i,:),A(i,:)] = prony(real(g_c2), 40, 40); toc
tic;    [B2(i,:),A2(i,:)] = stmcb(real(g_c2), 40, 40); toc
    % normalize filter to unity gain
    B(i,:) = B(i,:)./sum(B(i,:));
end
