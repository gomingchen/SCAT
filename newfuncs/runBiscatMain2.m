function sim = runBiscatMain2(varargin)
% RunBiscatMain  launches the BISCAT subprocesses using user specified options
%
% Simulation data is passed between subprocesses as a structure.
%
% SIM = RunBiscatMain runs BISCAT with the current settings in the variable
%       'cfg'.  If the variable is not found, it will attempt to load the
%       default settings in 'config.mat'.
% SIM = RunBiscatMain(CFGVAR) uses settings in the variable, CFGVAR.
% SIM = RunBiscatMain(CFGFILE) where CFGFILE is a string, uses settings in
%       the specified file, specified with a relative or absolute path name.
% SIM = RunBiscatMain(CFG,TIMESERIES) overrides the input file named in
%       'cfg.file_name' and runs BISCAT with TIMESERIES as the input to the
%       model.  TIMESERIES can be a filename or valid time series structure.
% SIM = RunBiscatMain(CFG,TIMESERIES,PANEL) executes up to the panel number
%       entered
%
% 'cfg' Structure Contents:
%
%                      GENERAL PARAMETERS
%           file_name: filename of time series input [string]
%      file_normalize: flag to normalize the signal energy [bool]
%      file_normlevel: desired RMS sound pressure level (db) [real]
%     file_timeseries: flag to display time series data [bool]
%       file_specgram: flag to display time series spectrogram [bool]
%    file_wignerville: flag to display time series pseudo WVD [bool]
%
%                      COCHLEAR BLOCK PARAMETERS
%          coch_panel: type of filterbank [string]
%            coch_ord: order of each filter [int]
%         coch_ripple: allowable ripple for some filter types [real]
%             coch_bw: -3dB bandwidth of each filter [real]
%        coch_fcenter: flag specifying constant Q/BW filterbank spacing [bool]
%           coch_fmin: minimum center frequency in filterbank [real]
%           coch_fmax: maximum center frequency in filterbank [real]
%          coch_steps: number of filters in bank [int]
%       coch_freqresp: flag to display frequecy response of filters [bool]
%            coch_bmm: flag to display basilar membrane motion plot [bool]
%    coch_cochleagram: flag to display cochleagram output [bool]
%
%                      NEURAL TRANSDUCTION PARAMETERS
%          neur_panel: type of neural response model [string]
%         neur_refper: refractory period used in model [real]
%         neur_reffcn: type of refractory function used in model [string]
%         neur_refval: value of relative refractory level, 0 to 1 [real]
%            neur_bmm: flag to display Basilar Membrane Motion plot [bool]
%
%          panel_rect: type of post-filtering rectification [string]
%            rect_bmm: flag to display Basilar Membrane Motion plot [bool]
%           panel_env: type of envelope smoothing filter [string]
%             env_ord: order of smoothing filter [int]
%              env_fc: -3dB cutoff frequency for LP smoothing filter [real]
%        env_freqresp: flag to display frequency response of filter [bool]
%
%                      ANALYSIS PARAMETERS
%
%
% 'ts' Structure Contents:
%
%                data: Nx1 or Nx2 array of time series signal [real]
%                time: Nx1 array of time vector (optional) [real]
%                  fs: constant sampling rate used in time series [real]
%           timestamp: date and time the input signal was saved [string]
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input error and parameter checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load functions in subdirectories into current path
% 
% subdirs = {'./sig_generate','./coch_filters', './neur_discharge', './neur_temporal', ...
%     './biscat_signals','./biscat_signals/hfm', './matlib_lib', ...
%     './delayline', './config_files', './gui_callbacks'};
% for str = subdirs
%     if isempty(strfind(path, str{:}))
%         addpath(str{:});
%     end
% end



switch nargin
    case 0
        % load default configuration settings
        warning('SCAT:InputArgs','No arguments passed, attempting to load default configuration file, config.mat')
        load('config.mat');
        
    case 1
        % load configuration settings
        if isstruct(varargin{1})
            cfg = varargin{1};
        elseif ischar(varargin{1}) && exist(varargin{1},'file')
            load(varargin{1});
        else
            error('Invalid input parameter used')
        end
        
    case 2
        % load configuration settings
        if isstruct(varargin{1})
            cfg = varargin{1};
        elseif ischar(varargin{1}) && exist(varargin{1},'file')
            load(varargin{1});
        else
            error('Invalid input parameter used')
        end
        
        % load time series data
        if isstruct(varargin{2})
            ts = varargin{2};
        elseif ischar(varargin{2}) && exist(varargin{2},'file')
            load(varargin{2});
        else
            error('Invalid input parameter used')
        end
        
    case 3
        % load configuration settings
        if isstruct(varargin{1})
            cfg = varargin{1};
        elseif ischar(varargin{1}) && exist(varargin{1},'file')
            load(varargin{1});
        else
            error('Invalid input parameter used')
        end
        
        % load time series data
        if isstruct(varargin{2})
            ts = varargin{2};
        elseif ischar(varargin{2}) && exist(varargin{2},'file')
            load(varargin{2});
        else
            error('Invalid input parameter used')
        end
        
        % mark panel number
        numPanel = varargin{3};
        
    otherwise
        error('Too many input arguments, not sure what to do!')
end

% fail if no configuration settings found
if ~exist('cfg','var')
    error('No configuration variable found named ''cfg''.')
end

% attempt to load data if not entered as parameter
if ~exist('ts','var')
    load(cfg.file_name)
end

% add time vector to data if not already present
if ~isfield(ts,'time')
    ts.time = (0:length(ts.data)-1)'./ts.fs;
end 

% ensure time series is a matrix of column vectors (Mx1 or Mx2)
if size(ts.data,1) < size(ts.data,2)
    ts.data = ts.data';
    %warning('RUNBISCATMAIN:DATA','Time series data should be stored as column vectors.  Transposing "ts.data"...');
end
if size(ts.time,1) < size(ts.time,2)
    ts.time = ts.time';
    %warning('RUNBISCATMAIN:DATA','Time series data should be stored as column vectors.  Transposing "ts.time"...');
end

% convert complex data to real
if ~isreal(ts.data)
    ts.data = real(ts.data);
end

% load('pulse.mat');
% 
% pe = [data'; ts.data(3501:end)]; % new pulse echo
% ts.data = pe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input signal(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize time series data
if cfg.file_normalize
    eRMS = sum(abs(ts.data).^2) / length(ts.data);
    ts.data = 10.^(cfg.file_normlevel/20) * ts.data / sqrt(eRMS);
end

% plot input signals
%plotTimeSeries(cfg,ts);

if ~exist('numPanel', 'var')
    numPanel = 3;
end

if (numPanel == 1); sim=[]; return; end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cochlear filter bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run cochlear filterbank process
sim.coch = runFiltBank(ts, cfg);

% plot cochlear filterbank signals
plotFiltBank_noplot(cfg,ts.fs,sim.coch);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neural Transduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run neural transduction process
% sim.neur = runNeurTran(sim.coch.bmm, ts, sim.coch.Fc, cfg);
% 
% % plot neural transduction signals
% plotNeurTran(cfg,ts,sim.coch.Fc,sim.neur);
% 
% if (numPanel == 2), return, end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cortical Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run neural analysis process
% sim.cn = runCochNuc(sim.neur, ts, sim.coch.Fc, cfg);
% 
% % plot neural analysis signals
% plotCochNuc(cfg,  sim.cn);
% name = varargin{2};
% ind = strfind(name, '_');
% ind2 = strfind(name, 'd');
% SNR = str2double(name(ind(1)+1:ind2-1));

% if (numPanel == 3), return, end

