
function [DF, d1] = categorizeEcho(echoL, K, Fs, Fc)
% after dechirping the spectrogram, categorize the dechirped echoes into
% separate echoes


target_M = K;
nchan = length(echoL);
Nind = length(K);
d1 = zeros(nchan, Nind);
timelen = 2.1; % how long a pulse is

    for p = 1:Nind
    gap = zeros(nchan,1);
    target_sm = target_M(p);
    for h = 1:nchan
        
        if isempty(echoL{h})
            gap(h) = NaN;
        else
            tem = echoL{h};
            kk = find((tem - target_sm)<timelen*Fs/1E3 & (tem - target_sm)>0);
            if isempty(kk)
                gap(h) = NaN;
            else
                if p > 1
                    if length(kk)>1
                        if tem(kk(1)) == d1(h,p-1)
                            gap(h) = tem(kk(2));
                        else
                            gap(h) = tem(kk(1));
                        end
                    else
                        if tem(kk(1)) == d1(h,p-1)
                            gap(h) = NaN;
                        else
                            gap(h) = tem(kk(1));
                        end
                        
                    end
                else
                    gap(h) = tem(kk(1));
                end
            end
        end
    end
    d1(:,p) = gap;
    end


    
    
% Plot dechirped image
    hold on
    Fc_backup = Fc;
%     scatter(zeros(size(Fc)),Fc);
    DF = zeros(nchan, Nind-1);
    for g = 2:4
        df = d1(:,g) - d1(:,1);
        kn = find(~isnan(df));
        mean_df = mean(df(kn));
        kkf = find(abs(df-mean_df)>500);
        if ~isempty(kkf)
            df(kkf) = NaN;
            Fc(kkf) = NaN;
        else
            Fc = Fc_backup;
        end
%         scatter(df, Fc);
        DF(:,g-1) = df;
    end
    
    
    
%     scatter(d1(:,3)-d1(:,1),Fc);
%     
%     scatter(d1(:,4)-d1(:,1),Fc);
%     scatter(d1(:,2)-d1(:,1),Fc);
%     scatter(d1(:,6)-d1(:,1),Fc);
%     scatter(d1(:,7)-d1(:,1),Fc);
%     scatter(d1(:,8)-d1(:,1),Fc);
end