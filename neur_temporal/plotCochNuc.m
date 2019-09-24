function plotCochNuc(cfg,cn)
% PLOTCOCHNUC  plots neural analysis signals
%
% PLOTCOCHNUC(CFG,CN)  is called from RunBiscatMain.m

% Plot synaptic transconducances (g_ex, g_in, g_re)
if cfg.cn_plottran
    figure;
    plot(cn.time, cn.g_ex, cn.time, cn.g_in)
    grid on;
    title(sprintf('Excitatory and Inhibitory Transconductances (N=%d, M=%d)', size(cn.spikes,2), cfg.cn_ntot))
    xlabel('Time (sec)')
    ylabel('Transconductance (dimensionless)')

    if cfg.cn_nrec > 0
        figure;
        plot(cn.time, cn.g_re)
        grid on;
        title(sprintf('Recurrent Synaptic Transconductances (N=%d, M=%d)', size(cn.spikes,2), cfg.cn_ntot))
        xlabel('Time (sec)')
        ylabel('Transconductance (dimensionless)')
    end
end

% Plot autocorrelation of output neurons
if cfg.cn_autocorr
    if sum(cn.spikes(:,1)) > 3
        [cn.A,cn.M] = autocorrspike(cn.V(:,1),cn.fs,0.2);
        figure;
        bar(cn.M*1e3, cn.A, 1);
        title('Autocovariance of single Bushy cell (AVCN) response')
        xlabel('Time (ms)')
        grid on
    else
        warning('Auto-Correlation Plot:  Not enough spikes present in CN output!')
    end
end

% Plot neuron membrane potential along with external stimulus vs. time
if cfg.cn_plotvmem
    figure;
    plot(cn.time, cn.V)
    grid on; hold on;
    plot([0 cn.time(end)],[cn.NN.Vth cn.NN.Vth],'--r')
    title(sprintf('Integrate and Fire Neuron Membrane Potential (N=%d, M=%d)', size(cn.spikes,2), cfg.cn_ntot))
    xlabel('Time (sec)')
    ylabel('Potential (V)')
    set(gca,'YLim',[cn.NN.Vreset*1.2 0.12])
end

% Plot ISI histogram of Bushy cell neurons
if cfg.cn_isihist
    if ~isempty(cn.ISI)
        figure;
        histpdf(cn.ISI,100,[0 max(cn.ISI)])
        title(sprintf('Cochlear Nucleus (CN) Interspike Interval (ISI) Histogram:  C_v = %.2f',cn.Cv))
    else
        warning('ISI Histogram:  No spikes present in CN output!')
    end
end

% display raster plot of neural activity
if cfg.cn_rasterplot
    figure;
    rasterplot(cn.spikes', cn.fs);
    title('Cochlear Nucleus Spike Generation')
end