function y = fb_genlatency(x,fs)
% 

% parameters
T = 500;            % integration time of energy estimate
K = -16e-6;         % amplitude latency delay constant
win = 'hamming';    % window function to use in energy estimate

% estimate energy in each channel
for i = 1:size(x,2)
    E(:,i) = energy(x(:,i), win, 1/T, T);
    %E(:,i) = db(E,'power');          % convert to dB
    %E(:,i) = E(:,i)./max(E(:,i));        % normalize to peak energy
end

if any(any(E == 0))
    E(E==0) = 1e-9;
end

% translate energy level to desired delay
deltau = (K*fs)*db(E,'power');       % group delay (converted to samples) -16us/dBA

%     % iteratively update filter weights and process data
%     N = size(B,2)-1;  % filter order
M=0;%M = size(A,2)-1;
L = length(x);
C = size(x,2);
D = ceil(max(max(deltau)));     % automatically size to largest desired delay

%     x = [zeros(D,size(x,2)); x];        % zero pad data vector
x = [zeros(C,D) x.'];
y = zeros(C,L);                 % preallocate memory for filter output
% 
%     %z = zeros(C,L);                   % preallocate memory for 2nd filter output
%     %e = zeros(C,L);                 % preallocate memory for energy estimate
%     %w = zeros(C,T);                   % preallocate memory for group delay filter
%     
%     % reverse coefficients for convolutions
%     %b = fliplr(B);
%     %a = fliplr(A);

for j=1:L
    % adjust filter coefficients (impulse response)
    w = sparse(C,D);
    w(sub2ind([C,D], (1:C), round(deltau(j,:)))) = 1;        % assign appropriate delay

    % perform group delay update on each channel
    y(:,j) = sum(w .* x(:,j:j+D-1),2);

end

y = y';

%     
%     % iterate over each time sample (all channels computed simultaneously)
%     for j=1:L
%         % "Direct Form II Transposed" implementation
%         y(:,j+M) = b * x(j:j+N);% - sum(a .* y(:,j:j+M),2);
%         
%         % Compute short-time energy estimate
%         % NEED TO CHANGE THIS TO Y, NOT X.. also, we're looking too far ahead!  Don't have future samples yet.
% %        e(1,j) = sum(win.*(x(j:j+T).^2));
% %        e(:,j) = sum(win.*(y(:,j:j+T).^2));
%         
%         % Adjust group delay filter coefficients
%         %w(:,j) = b_e./e(:,j);  % simple energy normalization MA filter - should see energy even out
%         
%         % Pass data through variable group delay filter
%         %z(:,j) = w(:,j:j+T) * x(j:j+T);
%         
%     end
%     coch.bmm = y(:,M+1:end).';       % remove initialization data and transpose
