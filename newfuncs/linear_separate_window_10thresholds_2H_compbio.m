
function [sd, sbr] = linear_separate_window_10thresholds_2H_compbio(wavParam)
% This function is to find the intercepts between one threshold value and
% pulse-echo pair.
% outputs: sd, the intercepts of echoes, without dechirping
%          sbr, the intercepts between thresholds with the broadcast,
%          without dechirping
%%% Retrieve wave parameters

sim = wavParam.simStruct;
NoT = wavParam.NoT; % which threshold, from 1 to 10
Fs = wavParam.Fs; % sampling frequency
twi = wavParam.callLenForMostFreq; % the length of echo/broadcast for most frequencies
twih = wavParam.callLenForHighFreq; % the length of echo/broadcast for high frequencies
twu = wavParam.callLenSpecial; % the minimum length of echo/broadcast for special frequencies, usually the lowest several frequencies
th_type = wavParam.th_type;
th_val = wavParam.startingThPercent;
sepbwBRand1stEcho = wavParam.SepbwBRand1stEchoinSmpls;
brStartSmpls = round(wavParam.whenBrStart.*Fs);
twb = wavParam.callLenBr;
c = wavParam.color;
coef_alt = wavParam.ALT;
NT = wavParam.NT; % number of threshold levels
% N1 = wavParam.twoHarmonicNNs(1);
% N2 = wavParam.twoHarmonicNNs(2);


stimIn = sim.coch.bmm;
N = size(stimIn, 2);

Fp = 10E3;
order = 30;
WLB = round(twb.*Fs);

lowpassFilt = lowpass_filt(Fs, Fp, order);

%%% For debugging purposes. change plotlog to 1 and plotlog2 to 0 in order
%%% to plot individual channels and the intercepts. Attention: that will
%%% produce many figures if left run for N times.

plotlog = 0;  % plot individual frequency channel, 1 for plot, 0 for not plot
plotlog2 = 1;  % plot the scatter plot for all frequency channels, 1 for plot, 0 for not plot
sbr = {};


% if NoT == 1 % for the first threshold, overwrite the input interval by setting higher/wider length of threshold.
%     twih = 1E-3;
%     twi = 1E-3;
% end

for i = 1:N
    
    %%% rectify the signal
    
    if i>60
        WL = twih*Fs;
    else
        WL = twi*Fs;
    end

    wf = stimIn(:,i);
    
    k = find(wf<0);
    wf(k) = 0;
    rect_wf = wf;
  
    tg = find(rect_wf<0);
    rect_wf(tg) = 0;
    sm_wf = lowpassFilt(rect_wf);
    hh = find(sm_wf<0);
    sm_wf(hh) = 0;
    [pmax, ~] = max(sm_wf(1:sepbwBRand1stEcho));
    SM_WF = (sm_wf)*100./(pmax);
    [newPmax, ~] = max(SM_WF);
    
    
    
    
    switch th_type
        case 'const'
            th_one = th_val;  % the lower boundary of thresholds across frequencies
        case 'varying' 
            th_one = 10/(i+th_val)*newPmax;  % the lower boundary of thresholds across frequencies
    end

    th_vec_pulse = linspace(th_one, 0.98*100, NT);
    thresh = th_vec_pulse(NoT);
    sm_wf = SM_WF(1:end);
    

    
    %%%%%%% FIGURE 1
    if plotlog
    figure, plot(sm_wf);
    frq_stp = (sim.coch.Fc(2) - sim.coch.Fc(1))/1E3; % in kHz
    title(sprintf('frequency %3.1f kHz', (sim.coch.Fc(1)/1E3 + (i-1)*frq_stp)));
    xlabel('time (ms)');
    ylabel('amplitude (norm.)');
    set(gca, 'XTick', 1:(5*Fs/1E3):length(sm_wf),'XTickLabel', 0:5:(length(sm_wf)/(Fs/1E3)),'FontName', 'Times', 'FontSize', 12);
    
    hold on,
    end
    %%%%%%% END FIGURE 1
    
%%% Find all values above threshold for the broadcast
    if i>20
        
        kk = brStartSmpls - 1 + find(sm_wf(brStartSmpls:sepbwBRand1stEcho-1)>=thresh);
        p_ind{i} = [];
        p_ind{i} = [p_ind{i} kk(1)];
            while ~isempty(kk)
                b = kk(1);
                kk = b + WLB - 1 + find(sm_wf(b+WLB:sepbwBRand1stEcho-1)>=thresh);
                if ~isempty(kk)
                    p_ind{i} = [p_ind{i} kk(1)]; % it's almost guaranteed to be not empty
                   
                end
            end
       
      gg = find(sm_wf(sepbwBRand1stEcho:end)>=thresh);
      e_ind = [];
      if ~isempty(gg)
          e_ind = sepbwBRand1stEcho + gg(1) - 1;
          a = sepbwBRand1stEcho+ gg(1) - 1;
      end
%%% Find the values above threshold for the echo     
            while ~isempty(gg)
                
                
                gg = a + WL - 1 + find(sm_wf(a+WL:end)>=thresh);
                if ~isempty(gg)
                    e_ind = [e_ind gg(1)]; % it's almost guaranteed to be not empty
                    a = gg(1);
                end
            end
            
            
    else % when i <= 20, only count the first intercept with the broadcast as broadcast
       WL2 = twu*Fs;
       e_ind = [];

       thresh = th_vec_pulse(NoT);
       kk = brStartSmpls - 1 + find(sm_wf(brStartSmpls:end)>=thresh);
       gg = find(sm_wf(sepbwBRand1stEcho:end)>=thresh);
       p_ind{i} = kk(1);
       
       
       if ~isempty(gg)
          e_ind = sepbwBRand1stEcho + gg(1) - 1;
          a = sepbwBRand1stEcho+ gg(1) - 1;
       end
       
       while ~isempty(gg)
               gg = a + WL2 - 1 + find(sm_wf(a+WL2:end)>=thresh);
                if ~isempty(gg)
                    e_ind = [e_ind gg(1)]; % it's almost guaranteed to be not empty
                    a = gg(1);
                end
       end 
        
    end
            
    % find the values in e_ind and p_ind that are too close to the
    % boundaries. 
    
    ke = find(abs(e_ind - sepbwBRand1stEcho) < 50);
    e_ind(ke) = [];
    
    GP = p_ind{i};
    
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
    
    
    %%% Account for amplitude latency trading effect
    
    % Calculate amplitude Latency Trading Effect
    if ~isfield(wavParam, 'ALToff') && ~isempty(e_ind)
        NE = length(e_ind);
    [time_in_us, loc] = amp_lat_trading_multiechoes(stimIn(:,i), 0.07, thresh, coef_alt, NE); % -: echo will be further delayed; + - 
    smpl_of_ALT = time_in_us.*1E-6*Fs;
    
    if smpl_of_ALT > 0
        smpl_of_ALT = 0;
    end
   
    
    if size(smpl_of_ALT,1)>1
        smpl_of_ALT = smpl_of_ALT';
    end
    kk = [];
    e_copy = e_ind;
    
    if ~isempty(e_copy) && ~isempty(smpl_of_ALT)
    for g = 1:length(smpl_of_ALT)
        gk = find(abs(e_copy - loc(g))<750);
        e_ind(gk) = e_ind(gk) + abs(floor(smpl_of_ALT(gk)));
    end
    end
    
    end
    

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
         
    else
        diff = [];
        
    end
    
    if ~isempty(e_ind)
        diff = e_ind;
        t_int = diff./Fs; % in second
        jj = find(t_int<twu);
        diff(jj) = [];
        sd{i} = diff;
    else
        sd{i} = [];
    end
    
    %%% PLOT 2 which produces the dechirped scatter image.
    if plotlog2

    hold on
    Fc = sim.coch.Fc;
    if ~isempty(diff)

        t_int = diff./Fs; % in second
        jj = find(t_int<twu);
        diff(jj) = [];
        
        %scatter(zeros(length(diff),1), Fc(i).*ones(1, length(diff)), ['o' c]);
        scatter(diff/Fs*1E3, Fc(i)*ones(1,length(diff)), ['o' c]); % choose the first pair to plot
%     else
%         scatter(zeros(1,1), Fc(i).*ones(1,1), ['o' c]);
    end
    
    if ~isempty(GP)
        sbr{i} = GP;
        scatter(GP/Fs*1E3, Fc(i).*ones(length(GP),1), 'om');
    else
        sbr{i} = [];
    end
    

    end

end

    
    
    
    