function neur = runNeurTran(bmm, ts, Fc, cfg)
% RUNNEURALDISCHARGE  Generates auditory spikes by various methods
%
% Currently available neural models:
%     'neur_rcf'        - Rectify, Compress, and Filter model
%     'neur_bio'        - Meddis (1986) IHC model
%     'neur_rnd'        - Random Poisson processes (ignores time series input)
%
% To be implemented:
%     'neur_hiel'       - Heil's onset spiking model as seen in Reijniers & Peremans
%


% Author:  Jason Gaudette
% Email:   jason.e.gaudette@navy.mil
% Date:    5/16/2011
% Version: 0.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run selected neural transduction model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.neur_panel
    case 'neur_rcf'
        neur = nt_rcf(ts.fs, bmm, Fc, cfg);     % Rectify, Compress, Filter Model
    case 'neur_bio'
        neur = nt_meddis(ts.fs, 1e4.*bmm, cfg);  % Biophysical Model
    case 'neur_rnd'
        neur = nt_random(ts, cfg);              % Random Spiking Process
    otherwise
        warning('NeurTran:MODE','Unknown neural discharge model:  %s', cfg.neur_panel)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike Statistical Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neur.ISI = genisi(neur.spikes(:),ts.fs);
neur.Cv = std(neur.ISI)/mean(neur.ISI);
