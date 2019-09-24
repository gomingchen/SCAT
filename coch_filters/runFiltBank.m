function coch = runFiltBank(ts, cfg)
% RUNFILTBANK  Runs time series data through a cochlear filterbank
%
% Currently available filterbank models:
%     'coch_butter'     - Butterworth filterbank
%     'coch_cheby1'     - Chebychev Type-I filterbank
%    *'coch_firls'      - FIR (least-squares) filterbank
%     'coch_gammatone'  - Gammatone filterbank
%    *'coch_gammachirp' - Gammachirp filterbank
%    *'coch_drnl'       - Dual-Resonant Non-Linear filterbank
%
% * = These functions not fully implemented
%

% Author:  Jason Gaudette
% Email:   jason.e.gaudette@navy.mil


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup center frequencies (Fc) and cutoffs (Fn,Wn) for filter bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.coch_fcenter
    case 1      % constant BW / linear period spacing
        ERB = cfg.coch_bw;
%         Tc = (0:1/(cfg.coch_steps-1):1)' .* (1/cfg.coch_fmax-1/cfg.coch_fmin) + 1/cfg.coch_fmin;
%         coch.Fc = 1./Tc;              % center frequency of filter
        stepFreq = (cfg.coch_fmax - cfg.coch_fmin)/(cfg.coch_steps-1);
        coch.Fc = (cfg.coch_fmin:stepFreq:cfg.coch_fmax);
        
    case 2      % constant Q
        % use default ERB method for now
        method = 'moore';
        switch lower(method)
            case{'lyon','stanley'}
                EarQ = 8;       % Lyon + Stanley Parameters (1988)
                minBW = 125;
                order = 2;
            case{'greenwood'}
                EarQ = 7.23824; % Greenwood Parameters (1990) as (nearly) in DSAM
                minBW = 22.8509;
                order = 1;
            case{'moore','glasberg'}
                EarQ = 9.26449; % Glasberg and Moore Parameters (1990)
                minBW = 24.7;
                order = 1;
            case{'wierddsam'}
                EarQ = 9.26; % Glasberg and Moore Parameters (1990)
                minBW = 15.719; %WORKS IF YOU SWCREW WITH THIS PARAMETER AS SO. . . 
                order = 1;
            otherwise
                error('Specified method not valid');
        end

        % compute ERB based on order and center frequency
        ERBlo = ((cfg.coch_fmin/EarQ)^cfg.coch_ord+ minBW^cfg.coch_ord) ^ (1/cfg.coch_ord);
        ERBhi = ((cfg.coch_fmax/EarQ)^cfg.coch_ord + minBW^cfg.coch_ord) ^ (1/cfg.coch_ord);
        overlap = (ERBhi/ERBlo)^(1/(cfg.coch_steps-1));
        ERB = ERBlo * (overlap.^(0:cfg.coch_steps-1))';
        coch.Fc = EarQ*(((ERB.^cfg.coch_ord) - (minBW.^cfg.coch_ord)).^(1/cfg.coch_ord));

    otherwise
        warning('SCAT:FiltBank','Invalid selection for filterbank center frequency spacing')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate cutoff frequency from center frequency, Fc, and bandwidth(s), ERB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coch.Fn = [coch.Fc-ERB/2 coch.Fc+ERB/2];
coch.Wn = 2*coch.Fn./ts.fs;         % normalize to Nyquist for digital filter

% create the filterbank using configuration settings
switch cfg.coch_panel
    case 'coch_butter'
        [B,A] = fb_butter(ceil(cfg.coch_ord/2), coch.Wn);
    case 'coch_cheby1'
        [B,A] = fb_cheby1(cfg.coch_ord, cfg.coch_ripple, coch.Wn);
    case 'coch_firls'
        B = fb_firls(cfg.coch_ord, coch.Wn);
        A = ones(cfg.coch_steps, 1);
    case 'coch_firwin'
        B = fb_firwin(cfg.coch_ord, coch.Wn);
        A = ones(cfg.coch_steps, 1);
    case 'coch_gammatone'
        [B,A] = fb_gammatone(ts.fs, coch.Fc, ERB);
    case 'coch_gammachirp'
        [B,A] = fb_gammachirp(ts.fs, coch.Fc, ERB, cfg.coch_fmax);
    case 'coch_drnl'
        error('SCAT:Cochlear_Block', 'DRNL filter not yet implemented!')
    otherwise
        error('SCAT:Cochlear_Block', sprintf('Unknown parameter: "%s"',cfg.coch_panel))
end

% save coefficients
coch.B = B;
coch.A = A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate basilar membrane movement for each cochlear frequency channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perform static filtering on each channel
coch.bmm = zeros(length(ts.data), cfg.coch_steps);        % preallocate memory
for i=1:cfg.coch_steps
    coch.bmm(:,i)=filter(B(i,:),A(i,:),ts.data);
end

% perform adaptive (time-variant) filtering for each time step
if cfg.coch_mode
    coch.bmm = fb_genlatency(coch.bmm, ts.fs);
end

coch.labels = (1 : floor(cfg.coch_steps/10) : cfg.coch_steps);    % setup default labels for BMM
