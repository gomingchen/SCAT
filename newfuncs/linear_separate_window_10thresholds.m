
function [sd, sgl,traces_echo] = linear_separate_window_10thresholds(wavParam)
% This function is to find the intercepts between one threshold value and
% pulse-echo pair.

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
c = wavParam.color;
coef_alt = wavParam.ALT;

stimIn = sim.coch.bmm;
N = size(stimIn, 2);

Fp = 10E3;
order = 30;
NT = 10; % number of threshold levels

lowpassFilt = lowpass_filt(Fs, Fp, order);

%%% For debugging purposes. change plotlog to 1 and plotlog2 to 0 in order
%%% to plot individual channels and the intercepts. Attention: that will
%%% produce many figures if left run for N times.

plotlog = 0;  % plot individual frequency channel, 1 for plot, 0 for not plot
plotlog2 = 1;  % plot the scatter plot for all frequency channels, 1 for plot, 0 for not plot
sgl = zeros(N,1);

traces_echo = zeros(N,1); % traces in the echo

if NoT == 1 % for the first threshold, overwrite the input interval by setting higher/wider length of threshold.
    twih = 1E-3;
    twi = 1E-3;
end

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
    MAX = max(rect_wf); % so you always have to ensure the broadcast is stronger than the waveform


    tt = find(rect_wf<0.01*MAX);
    rect_wf(tt) = 0.01*MAX;
    
    sm_wf = lowpassFilt(rect_wf);
    
    [pmax, ~] = max(sm_wf(1:sepbwBRand1stEcho));
    uu = find(sm_wf<0);
    sm_wf(uu) = 0.01*pmax;
    
    [pmin, ~] = min(sm_wf);
    SM_WF = (sm_wf - pmin)*100./(pmax-pmin);
    [newPmax, ~] = max(SM_WF);
    
    
    switch th_type
        case 'const'
            th_one = th_val;  % the lower boundary of thresholds across frequencies
        case 'varying' 
            th_one = 10/(i+th_val)*newPmax;  % the lower boundary of thresholds across frequencies
    end

    th_vec_pulse = linspace(th_one, 0.98*100, NT);
   
    sm_wf = SM_WF(1:end);
    
    %%% Calculate amplitude Latency Trading Effect
    time_in_us = amp_lat_trading(stimIn(:,i), 0.1, coef_alt); % -: echo will be further delayed; + - 
    smpl_of_ALT = time_in_us.*1E-6*Fs;
    
    if smpl_of_ALT > 0
        smpl_of_ALT = 0;
    end
    
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
        thresh = th_vec_pulse(NoT);
        kk = 50 - 1 + find(sm_wf(50:sepbwBRand1stEcho-1)>=thresh);
        p_ind{i} = [];
        p_ind{i} = [p_ind{i} kk(1)];
            while ~isempty(kk)
                b = kk(1);
                kk = b + WL - 1 + find(sm_wf(b+WL:sepbwBRand1stEcho-1)>=thresh);
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
            
            
    else
       WL2 = twu*Fs;
       e_ind = [];

       thresh = th_vec_pulse(NoT);
       kk = 50 - 1 + find(sm_wf(50:end)>=thresh);
        
       p_ind{i} = kk(1);
           while ~isempty(kk)
               a = kk(1);
               kk = a + WL2 - 1 + find(sm_wf(a+WL2:end)>=thresh);
               if ~isempty(kk)
                   e_ind = [e_ind kk(1)]; % it's almost guaranteed to be not empty
               end
           end 
        
    end
            
    % find the values in e_ind and p_ind that are too close to the
    % boundaries. 
    
    ke = find(abs(e_ind - sepbwBRand1stEcho) < 50);
    e_ind(ke) = [];
    
    GP = p_ind{i};
    kp = find(GP < 50);
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
    
    
    %%% Account for amplitude latency trading effect
    if ~isfield(wavParam, 'ALToff')
        e_ind = e_ind+abs(floor(smpl_of_ALT));
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
        
        if LE == 1
            if i < 6 
                traces_echo(i) = e_ind;
            else
                if i >= 43-25+1 && i <= 47-25+1
                    traces_echo(i) = e_ind;
                end
            end
        end
        % get a mean trace of later val(lower frequency) and earlier val(higher frequency)
        if i>47-25+1 && i<52-25+1 % in the region where the first and second harmonics overlap a bit
            kL = find(traces_echo(1:5));
            val_l = mean(traces_echo(kL));
            kH = find(traces_echo(19:23));
            val_h = mean(traces_echo(kH+19-1));
            if LE == 1 && LP == 2
                if (e_ind - val_l)>(e_ind - val_h)
                    % higher frequency, earlier values
                    sgl(i) = e_ind - GP(1);
                else
                    sgl(i) = e_ind - GP(2);
                    
                end
            else
                sgl(i) = e_ind(1) - GP(1);
            end
            
        else
         
            sgl(i) = e_ind(1) - GP(1);
        end
         
    else
        diff = [];
        sgl(i) = NaN;
        
    end
    
    if ~isempty(e_ind)
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
        
        scatter(zeros(length(diff),1), Fc(i).*ones(1, length(diff)), ['o' c]);
        scatter(diff/Fs*1E3, Fc(i)*ones(1,length(diff)), ['o' c]); % choose the first pair to plot
    else
        scatter(zeros(1,1), Fc(i).*ones(1,1), ['o' c]);
    end
    end
    

end

end

    
    
    
    
