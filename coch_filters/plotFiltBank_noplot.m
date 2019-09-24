function plotFiltBank_noplot(cfg,fs,coch)
% PLOTFILTBANK  plots cochlear filterbank signals
%
% PLOTFILTBANK(CFG,COCH)  is called from RunBiscatMain.m

numCh = size(coch.B,1);


% plot filter frequency responses
pDelay = 1e-3;    % plot delay between frequency channels
% if cfg.coch_tfplot
%     
%     % plot resulting filter output per channel
%     figure;
%     cmap = jet(numCh);      % define colormap
% 
%     for i=1:numCh
%         
%         % plot magnitude response
%         [H,F]=freqz(coch.B(i,:), coch.A(i,:), 4096,fs);
%         subplot(2,1,1);
%         plot(F./1000, db(abs(H)), 'color', cmap(i,:));
%         axis([0 fs/2000 -80 3]);
%         grid on; hold on;
%         title(sprintf('f_c=%.1f, BW=%.1f',coch.Fc(i),diff(coch.Fn(i,:))));
%         xlabel('Frequency (kHz)');
%         ylabel('Magnitude (dB)');
% %        title(cfg.panel_filt);
% 
%         % plot phase response
%         subplot(2,1,2);
%         plot(F./1000, 180.*unwrap(angle(H))./pi, 'color', cmap(i,:));
% %        axis([0 fs/2000 0 275])
%         grid on; hold on;
%         xlabel('Frequency (kHz)')
%         ylabel('Phase (degrees)')
% 
%         pause(pDelay)
%     end
% end

% plot filter group delay response
% if cfg.coch_grpdelay
%     figure;
%     cmap = jet(numCh);      % define colormap
% 
%     npts = 4096;
%     F = linspace(0,fs/2,npts+1);
%     F = F(1:end-1);
%     
%     % plot group delay response of each channel
%     maxGD = 0;
%     for i=1:numCh
%         GD = grpdelay(coch.B(i,:), coch.A(i,:), npts, fs);
%         plot(F./1000, GD, 'color', cmap(i,:));
%         grid on; hold on;
%         maxGD = max(maxGD,max(GD(50:end-50)));
%         axis([0 fs/2000, 0 maxGD*1.1]);
%         pause(pDelay)
%     end
%     title('Filter Group Delay Response')
%     xlabel('Frequency (kHz)')
%     ylabel('Group Delay (Samples)')
% end

% plot filter phase delay response
%%%%%%%%% NEED TO ADD THIS TO GUI %%%%%%%%%%
if 0 %cfg.coch_phsdelay
    figure;
    cmap = jet(numCh);      % define colormap

    npts = 4096;
    F = linspace(0,fs/2,npts+1);
    F = F(1:end-1);
    
    % plot phase delay response of each channel
    %maxPD = 0;
    for i=1:numCh
        [PD,w] = phasedelay(coch.B(i,:), coch.A(i,:), npts, fs);
        plot(F./1000, PD, 'color', cmap(i,:));
        grid on; hold on;
        %maxPD = max(maxPD,max(PD(50:end-50)));
        %axis([0 fs/2000, 0 maxPD*1.1]);
        pause(pDelay)
    end
    title('Filter Phase Delay Response');
    xlabel('Frequency (kHz)');
    ylabel('Phase Delay (degrees)');
end


% plot basilar membrane movement for each frequency channel
if cfg.coch_bmmplot
%     figure;
    BMMplot(coch.bmm, fs, coch.Fc, coch.labels);
    title(sprintf('%s', cfg.coch_panel),'Interpreter','none');
end
