
function [sd, sgl] = linear_separate_window_10thresholds_2H_shizuko(wavParam)
%%% To find the intercepts between the threshold and the smoothed waveform
%%% at each frequency channel after the broadcast-echo pair goes through
%%% the filter bank

%%% INPUT: wavParam, the structure that contains the parameters of general
%%% parameters/descriptions of the bmm
%%% OUTPUT: sd, cell array that contains the intercepts of broadcast and
%%% all echoes
%%%         sgl, array that contains only the intercepts of the first echo


%%% Assign the parameters
sim = wavParam.simStruct;
NoT = wavParam.NoT; % which threshold, from 1 to 10
Fs = wavParam.Fs; % sampling frequency
twi = wavParam.callLenForMostFreq; % the length of echo/broadcast for most frequencies
twih = wavParam.callLenForHighFreq; % the length of echo/broadcast for high frequencies
twu = wavParam.callLenSpecial; % the minimum length of echo/broadcast for special frequencies, usually the lowest several frequencies
br_startSec = wavParam.whenBrStart;
th_type = wavParam.th_type;
th_val = wavParam.startingThPercent;
sepbwBRand1stEcho = wavParam.SepbwBRand1stEcho;
c = wavParam.color;
stimIn = sim.coch.bmm;
N = size(stimIn, 2);
sepFlag = wavParam.sepFlag;

Fp = 10E3;
order = 30;
NT = 10; % number of threshold levels

lowpassFilt = lowpass_filt(Fs, Fp, order);

plotlog = 0;  % plot individual frequency channel, 1 for plot, 0 for not plot
plotlog2 = 1;  % plot the scatter plot for all frequency channels, 1 for plot, 0 for not plot
sgl = zeros(N,1);

for i = 1:N
    %%% rectify the signal
    if i>60
        WL = round(twih*Fs);
    else
        WL = round(twi*Fs);
    end

    wf = stimIn(:,i);
    
    k = find(wf<0);
    wf(k) = 0;
    rect_wf = wf;
    MAX = max(rect_wf);

    tt = find(rect_wf<0.01*MAX);
    rect_wf(tt) = 0.01*MAX;
    
    sm_wf = lowpassFilt(rect_wf);
    [pmax, ~] = max(sm_wf);
    uu = find(sm_wf<0);
    sm_wf(uu) = 0.01*pmax;
    
    [pmin, ~] = min(sm_wf);
    SM_WF = (sm_wf - pmin)*100./(pmax-pmin);
    [newPmax, ~] = max(SM_WF);
    
    brStartSam = round(br_startSec*Fs);
    
    switch th_type
        case 'const'
            th_one = th_val;  % the lower boundary of thresholds across frequencies
        case 'varying' 
            th_one = 10/(i+th_val)*newPmax;  % the lower boundary of thresholds across frequencies
    end

    th_vec_pulse = linspace(th_one, 0.98*100, NT);
    sm_wf = SM_WF(1:end);
    
%%% find the sample number that separates the broadcast and the following echoes
    if sepFlag
        Sep = wavParam.Sep.*Fs;
    else
        [Sep, yu] = findSep(sm_wf, 2, order); % it's second peak
    end

    
    %%%%%%% FIGURE 1
    if plotlog
    figure, plot(sm_wf);
    title(sprintf('frequency %d kHz', (sim.coch.Fc(1)/1E3 + i - 1)));
    xlabel('time (ms)');
    ylabel('amplitude (norm.)');
    set(gca, 'XTick', 1:2*Fs/1E3:length(sm_wf),'XTickLabel', 0:2:(length(sm_wf)/(Fs/1E3)),'FontName', 'Times', 'FontSize', 12);
    ylim([0,150]);
        if ~sepFlag
            hold on,
            plot(yu,'r');
            scatter(Sep, yu(Sep), 'or');
        end
    end
    %%%%%%% END FIGURE 1

    %%%%%%% FIND THE INTERCEPTS OF BROADCAST AND ECHOES
    %%% BROADCAST
      thresh = th_vec_pulse(NoT);
      kk = brStartSam - 1+find(sm_wf(brStartSam:end)>=thresh);
      p_ind{i} = [];
      if ~isempty(kk)
          p_ind{i} = [p_ind{i} kk(1)];
      end

    %%% ECHOES
      gg = find(sm_wf(Sep:end)>=thresh);
      e_ind = [];
      if ~isempty(gg)
          e_ind = Sep + gg(1) - 1;
          a = Sep+ gg(1) - 1;
      end
  
      while ~isempty(gg)
                
           gg = a + WL - 1 + find(sm_wf(a+WL:end)>=thresh);
           if ~isempty(gg)
                e_ind = [e_ind gg(1)]; % it's almost guaranteed to be not empty
                a = gg(1);
           end
      end

            

    GP = p_ind{i};
    ke = find(abs(e_ind - GP(1))<sepbwBRand1stEcho*Fs/1E3);
    e_ind(ke) = [];
    kp = find(GP < br_startSec*Fs);
    GP(kp) = [];
    LP = length(GP);
    LE = length(e_ind);
    
    %%%% FIGURE 2
    if plotlog
    hold on
    scatter(GP, sm_wf(GP), 'filled');
    Q = e_ind;
    scatter(Q, sm_wf(Q), 'filled');
    end
    %%%% END FIGURE 2
    

    if LE ~= 0
        
        if LP == 1 || LP == LE
            diff = e_ind - GP;
        else
            diff = zeros(LP, LE);
            for j = 1:LP
                diff(j,:) = e_ind - GP(j);
            end
            diff = diff(:);
        end
        
    
        if ~isempty(GP)
            sgl(i) = e_ind(1) - GP(1);
        else
            sgl(i) = NaN;
        end
         
    else
        diff = [];
        sgl(i) = NaN;
        
    end
    

    sd{i} = [GP e_ind];
    
    %% the plot on the dechirped spectrogram
    if plotlog2

    hold on
    Fc = 1E-3.*sim.coch.Fc;
    ylabel('frequency (kHz)');
    xlabel('time (ms)');
    
    if ~isempty(diff)
        diff = e_ind;
        t_int = diff./Fs; % in second
        jj = find(t_int<twu);
        diff(jj) = [];
        if ~isempty(GP)
            scatter(GP(1)/Fs*1E3, Fc(i), 'or');
        end
        scatter(diff/Fs*1E3, Fc(i)*ones(1,length(diff)), ['o' c]); % choose the first pair to plot
    else
        scatter(e_ind/Fs*1E3, Fc(i).*ones(length(e_ind),1), ['o' c]);
    end
    end
    
     
end

end

    
    
    
    