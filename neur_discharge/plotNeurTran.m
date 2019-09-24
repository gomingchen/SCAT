function plotNeurTran(cfg,ts,Fc,neur)
% PLOTNEURTRAN  plots neural transduction signals
%
% PLOTNEURTRAN(CFG,TS,NEUR)  is called from RunBiscatMain.m

% plot panel specific signals
switch cfg.neur_panel
    case 'neur_rcf'
        % plot smoothing filter frequency response
        if cfg.neur_rcf_freqresp && cfg.neur_rcf_filt ~= 2
            F=10.^(3:.02:log10(ts.fs/2));
            H=freqz(neur.filt_b, neur.filt_a, F, ts.fs);
            figure; subplot(2,1,1);
            semilogx(F,db(abs(H)));
            title(sprintf('Smoothing filter frequency response (%gkHz, N=%d)', cfg.neur_rcf_fc./1000, cfg.neur_rcf_ord),'Interpreter','none');
            xlabel('Frequency (Hz)')
            ylabel('Magnitude (dB)')
            axis([F(1) F(end) -90 10])
            grid on;
            subplot(2,1,2);
            semilogx(F,angle(H)*180/pi);
            set(gca,'XLim',[F(1) F(end)])
            ylabel('Phase (degrees)')
            xlabel('Frequency (Hz)')
            grid on;
        end
        
    case 'neur_bio'
        % plot internal biophysical states from Meddis model
        if cfg.neur_biol_plotstates
            figure; mesh(1e3*ts.time, (1:cfg.coch_steps), neur.k');   %1e-3*sim.coch.Fc
            xlabel('Time (ms)');
            ylabel('Frequency Channel');
            title('k(t) = membrane permeability')
            colorbar

            figure; mesh(1e3*ts.time, (1:cfg.coch_steps), neur.c');
            xlabel('Time (ms)');
            ylabel('Frequency Channel');
            title('c(t) = neurotransmitter in cleft')
            colorbar

            figure; mesh(1e3*ts.time, (1:cfg.coch_steps), neur.q');
            xlabel('Time (ms)');
            ylabel('Frequency Channel');
            title('q(t) = neurotransmitter pool')
            colorbar

            figure; mesh(1e3*ts.time, (1:cfg.coch_steps), neur.w');
            xlabel('Time (ms)');
            ylabel('Frequency Channel');
            title('w(t) = reprocess store')
            colorbar
        end
end


% plot probability of event (P_spike)
if cfg.neur_probplot
    figure;
%    imagesc(1e3*ts.time, 1:size(neur.prob,2), neur.prob');
    plot(1e3*ts.time, neur.prob');
    xlabel('Time (ms)')
    ylabel('Frequency Channel')
    title('P_{Event}')
    set(gca,'ydir','normal')
    map = colormap;
    map(1,:) = [0 0 0];
    colormap(map);
    caxis([0 .025])

end


% plot interspike-interval (ISI) histograms
if cfg.neur_isihist
    if ~isempty(neur.ISI)
        figure; histpdf(neur.ISI,100,[0 max(neur.ISI)]);
        title(sprintf('Auditory Nerve Cell (ANC) Interspike Interval (ISI) Histogram:  C_v = %.2f',neur.Cv))
        xlabel('Time (sec)')
        ylabel('ISI Frequency')
    else
        warning('PLOT:NoData','ISI Histogram:  No spikes present in Neural Transduction output!')
    end
end


% display raster plot of neural activity
if cfg.neur_rasterplot
    figure;
    rasterplot(1e3*ts.time,Fc*1e-3,neur.spikes');
    title('Auditory Nerve Cell Spike Generation');
end

