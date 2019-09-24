function [V,g_ex,g_in,g_re] = leakyifsyn(T, N, S, W, R, Fs)
% LEAKYIFSYN  Simulates network of N homogenous integrate and fire neurons
%
% The membrane potential of each neuron is simulated using a leaky
% integrate and fire neural model.  The model updates membrane potential
% and synaptic conductances at each increment in time (a time-based
% simulation, not event-based).  Weights of synapses are fixed and are set
% by the supplied synapse interconnection matrix. Synaptic input is
% incremented postsynaptically using exponential or alpha functions.
%
% [V,Gex,Gin,Gre] = LEAKYIFSYN(T, N, S, W, Fs) returns a MxN matrix
%   of the time-dependent membrane potentials for each neuron, N.
%   Input:
%     T is an mx1 vector of simulation time
%     N is a struct containing the following named neuron parameters:
%         Vreset, Vrest, Vth, tau_m, tau_ex, tau_in, tau_re
%     S is an mxn matrix of presynaptic applied current
%     W is the synaptic interconnection matrix.  This pxn matrix contains
%       weights in the closed interval [-1 1].  p is the number of IF
%       neurons and n is the number of presynaptic inputs.
%     R is a recurrent synaptic interconnection matrix.  This pxp matrix
%       contains weights as above, but to/from each IF neuron.  Each row,
%       i, describes the synaptic weight from neuron i to the other p-1
%       neurons.  I.e. R(2,4) = 1 means neuron 2 propagates a spike to
%       neuron 4 with weight 1.  The main diagonal should always be zeros!
%     Fs is a scalar of the sampling rate

dt = 1/Fs;                  % calculate time increment

% initialize output matrices
V = N.Vrest*ones(length(T),size(W,1));  % initialize V to initial condition
g_ex = zeros(length(T),size(W,1));      % initialize g_ex to 0
g_in = zeros(length(T),size(W,1));      % initialize g_in to 0
g_re = zeros(length(T),size(W,1));      % initialize g_re to 0

% internal neuron model parameters
%g_exDel = 0.015;            % define incremental value of g_ex after each presynaptic AP
%g_inDel = 0.05;             % define incremental value of g_in after each presynaptic AP
%g_reDel = 0.05;              % define incremental value of g_re after each recurrent AP
E_ex = 0;                   % excitatory reversal potential of 0mV
E_in = -0.07;               % inhibitory reversal potential of -70mV

% Calculate the membrane potential over time iteratively
for t=3:length(T)
    
    % Iterate over each neuron (rows of matrix W)
    for n=1:size(W,1)
        
        % update membrane potential with derivative of voltage directly (5.8)
        if V(t-1,n) < N.Vth
            
            % implement differential equations as described in Song, Miller, and Abbott (2000)
            V_ex = g_ex(t-1,n)*(E_ex - V(t-1,n));  % excitatory synaptic input (modeled as conductance)
            V_in = g_in(t-1,n)*(E_in - V(t-1,n));  % inhibitory synaptic input (modeled as conductance)
            V_re = g_re(t-1,n)*(E_ex - V(t-1,n));  % recurrent synaptic input (modeled as conductance)
            dVdt = (N.Vrest - V(t-1,n) + V_ex + V_in + V_re) / N.tau_m;     % sum each voltage contribution and multiply by membrane time constant (5.8)
            V(t,n) = V(t-1,n) + dt*dVdt;              % incrementally add dV to last V

            % increment external synaptic conductances using exponential decay to zero
            g_ex(t,n) = g_ex(t-1,n) * (1 - dt/N.tau_ex);
            g_in(t,n) = g_in(t-1,n) * (1 - dt/N.tau_in);
            
            % increment recurrent synaptic conductances using exponential decay
            g_re(t,n) = g_re(t-1,n) * (1 - dt/N.tau_re);

        else
            % generate spike and reset potential
            V(t-1,n) = 0.1;               % set last value to spike
            V(t,n) = N.Vreset;            % set current value to Vreset

        end

        % detect presynaptic spikes and update conductances
        sEvent = S(t-1,:) .* W(n,:);
        if any(sEvent)

            % increase g_ex/in by predefined amount
            g_ex(t,n) = g_ex(t,n) + N.g_exDel * sum(sEvent(sEvent>0));
            g_in(t,n) = g_in(t,n) + N.g_inDel * sum(sEvent(sEvent<0));
            
        end
        
        % detect recurrent spikes and update conductances
        rEvent = [V(t-2,:)>0] .* R(:,n)';
        if any(rEvent)
            
            % propagate spike through to recurrent network using weights R
            g_re(t,n) = g_re(t,n) + N.g_reDel * sum(R(n,find(rEvent)));
            
        end
        
    end
end