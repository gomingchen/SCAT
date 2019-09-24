function temp = runTemporal(neur, ts, Fc, cfg)
% RUNTEMPORAL  Processes the auditory nerve cell (ANC) output in the
% temporal domain by the selected method
%
% Currently available neural models:
%     'temp_cn' - Cochlear nucleus (Integrate-and-fire) neural model
%     'neur_sc' - Spectrogram Correlation (SC) model from SCAT
%


% Author:  Jason Gaudette
% Email:   jason.e.gaudette@navy.mil
% Date:    6/7/2011
% Version: 0.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run selected temporal processing model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.neur_panel
    case 'neur_cn'
        temp = temp_cn(ts.fs, neur, Fc, cfg);   % Cochlear Nucleus model
    case 'neur_sc'
        temp = temp_sc(ts.fs, neur, Fc, cfg);   % Spectrogram Correlation model
    otherwise
        warning('NEURTEMP:MODE','Unknown temporal processing model:  %s', cfg.neurtemp_panel)
end
