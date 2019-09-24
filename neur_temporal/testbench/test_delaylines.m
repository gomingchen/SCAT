% construct delay lines comparing both, standard and sparse matrices

clc
clear

% parameters
fs = 1e6;       % "sampling" rate of delay lines
dMax = 0.01;    % delay line length
tMax = 0.1;     % simulation length
mode = 1;

% matrix dimensions
M = 100;                % number of CN neurons
N = ceil(dMax*fs);      % number of NLL neurons (per frequency channel)


% initialize delay lines
D0 = sparse(M,N,false);         % sparse matrix
D1 = false(M,N);                % logical matrix


% initial event
D0(:,1) = true; % D0(:,1:M) = speye(M);
D1(:,1) = true; % D1(:,1:M) = eye(M);


% define image axes
T = 1:N;
F = 1:M;

% update delay lines to propagate spike events
%figure(1)


tic
switch mode

    %% shift by indexing
    case 1
        
        % 5.5 sec on Dell for fs = 1e6, tmax = 0.1, dmax = 0.01, M = 100
        % 22 sec on MBP
        I0 = sparse(1:M,1,false);  I0(50) = true;
        for i = 1:round(tMax*fs)
            D0(:,i) = I0;
            image(T,F,D0)
        end
        
    case 2
        
        % > 5 min on Dell for fs = 1e6, tmax = 0.1, dmax = 0.01, M = 100
        % 1117 sec (18.6 min) on MBP
        I1 = false(M,1);
        for i = 1:round(tMax*fs)
            D1(:,i) = I1;
            %image(T,F,D1)
        end
    
    %% shift by matrix replication
    case 3
        % 169 sec on Dell for fs = 1e6, tmax = 0.1, dmax = 0.01, M = 100
        % 35 sec on MBP
        %image(T,F,D0);
        %colormap([1 1 1;0 0 0])
        I0 = sparse(1:M,1,false);
        for i = 1:round(tMax*fs)
            D0 = [I0 D0(:,1:end-1)];
        %    image(T,F,D0);
        %    pause(0.001)
        end
        
    case 4
        % 447 sec on Dell for fs = 1e6, tmax = 0.1, dmax = 0.01, M = 100
        % 414 sec on MBP
        %image(T,F,D0);
        %colormap([1 1 1;0 0 0])
        I1 = false(M,1);
        for i = 1:round(tMax*fs)
            D1 = [I1 D1(:,1:end-1)];
        %    image(T,F,D1);
        %    pause(0.001)
        end
        
    %% multiplication by NxN shift matrix
    case 5
        % 25 sec on MBP
        S0 = speye(N);
        S0 = [zeros(N,1) S0(:,1:end-1)];
        for i = 1:round(tMax*fs)
            D0 = D0 * S0;
        end

        
    case 6
        % out of memory on MBP
        S1 = eye(N);
        S1 = [zeros(N,1) S1(:,1:end-1)];
        for i = 1:round(tMax*fs)
            D1 = D1 * S1;
        end

end

toc
