function [p,c,q,w,k] = ANC_Meddis(S,fs,varargin)
% ANC_Meddis  Implementation of Ray Meddis' (1991) inner hair cell model
%
% p = ANC_Meddis91(S,fs) - uses all default model parameters
% [p,c,q,w,k] = ANC_Meddis91(S,fs,A,B,g,y,l,r,x) - allows more flexibility
% 
% Input:
%   S - signal(s) to be evaluated
%   fs - sampling rate of signals
%   A,B,g,y,l,r,x - IHC model parameters (refer to paper for more details)
%     A - offset for membrane permeability, k, to induce random spikes
%     B - adjusts curvature of membrane permeability, k, to input stimuli
%     g - gain of overall membrane permeability, k
%     y - replenish rate of NT in pool, q, by factory
%     l - loss rate of NT in cleft, c
%     r - reuptake rate of NT in cleft, c
%     x - reprocess rate of NT in pool, q, from store, w
%
% Output:
%   p - Probability of a spike occuring
%   c - Quantity of neurotransmitter in synaptic cleft
%   q - Free transmitter pool
%   w - Reprocessing store
%   k - Membrane permeability

% Note:  Implementation is based on HUTear2 mex functions

% ensure signal is column vector
if size(S,2) == length(S)
    S=S';
end

% define default model parameters
m=1;        % maximum capacity of neurotransmitter in pool, q(t)
h=500000;   % multiplication factor to produce spike probabality, p(t)

A=2;
B=5000;
g=200;
y=1.05;
l=3500;
r=6580;
x=66.31;

% parse optional parameters
if nargin == 3
    opt = varargin{1};
    A=opt(1);
    B=opt(2);
    g=opt(3);
    y=opt(4);
    l=opt(5);
    r=opt(6);
    x=opt(7);
end

% assign convenience factors
dt=1/fs;
gdt=g*dt;
ydt=y*dt;
ldt=l*dt;
rdt=r*dt;
xdt=x*dt;
hdt=h*dt;

% initialize state variables
k = zeros(size(S));
c = zeros(size(S));
q = zeros(size(S));
w = zeros(size(S));
p = zeros(size(S));

% % set initial conditions
% k(1,:)=g*A/(A+B);
% c(1,:)=k(1,:)*y*m/(y*(l+r)+k(1,:)*l);
q(1,:) = m;  %c(1)*(l+r)/k(1);
% w(1,:)=c(1)*r/x;

% step through time incrementally
for t=2:length(S)

    % iterate over each IHC/frequency channel
    for n=1:size(S,2)

        % calculate nonlinear BMM based on instantaneous amplitude of signal
        k(t,n) = gdt * max(S(t-1,n)+A,0) / (S(t-1,n)+A+B);

        % calculate model variables
        replenish = max(ydt*(m-q(t-1,n)), 0);
        eject = k(t-1,n)*q(t-1,n);
        loss = ldt*c(t-1,n);
        reuptake = rdt*c(t-1,n);
        reprocess = xdt*w(t-1,n);

        % implement differential equations of neurotransmitter model
        q(t,n) = q(t-1,n) + replenish - eject + reprocess;
        c(t,n) = c(t-1,n) + eject - loss - reuptake;
        w(t,n) = w(t-1,n) + reuptake - reprocess;

        % calculate the running probability of action potential
        p(t,n) = hdt*c(t,n);
    end
end
