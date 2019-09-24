function neur = nt_meddis(fs, bmm, cfg)
% NT_MEDDIS  Runs Meddis' 1986 IHC model on the basilar membrane motion
%
% neur = NT_MEDDIS(FS,BMM,CFG)

%%% TODO:  MAKE USE OF SPARSE ARRAYS TO SAVE MEMORY!
    
%% ----------- Bugs -----------
% R only calculated for HSR, needs to be included for MSR & LSR processes
% Returns neur.hsr, should remove temporary values and only include
% neur.prob & neur.spikes


% assign input parameters
Meddis_params = numeric(cfg.neur_biol_parm);
HSR_fanout = Meddis_params(1,1);        %numeric(cfg.neur_biol_hsrparam);
MSR_fanout = Meddis_params(2,1);        %numeric(cfg.neur_biol_msrparam);
LSR_fanout = Meddis_params(3,1);        %numeric(cfg.neur_biol_lsrparam);
HSR_params = Meddis_params(1,2:end);    %numeric(cfg.neur_biol_hsrfanout);
MSR_params = Meddis_params(2,2:end);    %numeric(cfg.neur_biol_msrfanout);
LSR_params = Meddis_params(3,2:end);    %numeric(cfg.neur_biol_lsrfanout);
aPer = numeric(cfg.neur_biol_absper)*1e-3;
rPer = numeric(cfg.neur_biol_relper)*1e-3;
c_r = numeric(cfg.neur_biol_cr);


% calculate final matrix dimensions
T = size(bmm,1);        % get time dimension
N = size(bmm,2);        % get channel dimension
L = LSR_fanout + MSR_fanout + HSR_fanout;       % total fanout
Ntot = L * N;     % number of column vectors for spike matrix


% generate ANC output from Meddis (1991) model (output is P[spike]/dt)
if HSR_fanout
    % init spike event matrix
    neur.hsr.spikes = sparse(T, N*HSR_fanout, false);
    %neur.hsr.spikes = boolean(zeros(T,N*HSR_fanout));
    
    if cfg.neur_biol_plotstates
        % init internal state vectors
        neur.k = zeros(T,N);
        neur.c = zeros(T,N);
        neur.q = zeros(T,N);
        neur.w = zeros(T,N);
        
        [neur.hsr.prob, neur.c, neur.q, neur.w, neur.k] = ...
            ANC_Meddis(bmm, fs, HSR_params);
    else
        neur.hsr.prob = ANC_Meddis(bmm, fs, HSR_params);
    end
    
    % expand probabilities with fanout
    neur.hsr.prob = kron(neur.hsr.prob, ones(1,HSR_fanout));        % UPDATE TO USE KRONPROD 3RDPARTY FUNCTION - THIS IS TOO SLOW
    
    % init uniform random process
    R = rand(T,N*HSR_fanout);

end

if MSR_fanout
    % init spike event matrix
    neur.msr.spikes = sparse(T, N*MSR_fanout, false);
    
    neur.msr.prob = ANC_Meddis(bmm, fs, MSR_params);
end

if LSR_fanout
    % init spike event matrix
    neur.lsr.spikes = sparse(T, N*LSR_fanout, false);
    
    neur.lsr.prob = ANC_Meddis(bmm, fs, LSR_params);
end


% generate spikes based on P[spike] at each time instant
if (aPer > 0) || (rPer > 0)
    for n = 1:Ntot     % iterate over each column vector in Meddis probability matrix
        idx = 1;
        % iterate over each spike in time
        while true
            idx = idx + find(neur.hsr.prob(idx:end,n) > R(idx:end,n), 1) - 1;
            if isempty(idx), break, end

            % mark spike event and update P[spike] for refractoriness
            neur.hsr.spikes(idx,n) = 1;

            % implement meddis91 refractory effects of transmitter available
            neur.hsr.prob(idx:end,n) = Ref_Meddis(neur.hsr.prob(idx:end,n), fs, aPer, rPer, c_r);
        end
    end
else
    neur.hsr.spikes(neur.hsr.prob > R) = 1;
end


% TEMPORARY FIX
neur.spikes = neur.hsr.spikes;
neur.prob = neur.hsr.prob;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = numeric(param)
% converts string input parameters to numeric type

if ischar(param)
    if numel(param) > 1
        param = str2num(param); %#ok<ST2NM>
    else
        param = str2double(param);
    end
end