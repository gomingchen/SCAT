
% for two harmonics, there are overlaps in frequency range between 40-60
% kHz.
% inputs, br, the broadcast cell array
%         echosq, the echo sequence cell array
%         R1, the reference (E1) for E1,
%         R2, the reference (E2) for E2
% outputs, D1, the dechirped sequences for first harmonic
%          D2, the dechirped sequence for 2nd harmonic
function [D1, D2, b] = dechirping(brdcst, echoL, R1, R2, Necho, echoLen, plotflag)


%load('testmat_dechirp');
%Necho = 10;
N = 10; % only count 25kHz and up
Nc = length(brdcst);
% for broadcast
u = 0;
b = zeros(Nc-N,2);
hb = 1; % harmonics indicator
he = 1; % harmonics flip for echoes
E = zeros(Nc-N, Necho, 2);

flipB = 0;
for x = N+1:Nc
    if length(brdcst{1,x}) == 2
        flipB = 1;
        break;
    end
end

% flipB is the indicator for if length(brdcst{1,anycell}) == 1

if flipB == 0
    u = 0;
    for a = N+1: Nc
        u = u+1;
        d = brdcst{1,a} - brdcst{1,a-1};
        if abs(d) > echoLen
            hb = 2;
        end
        b(u,hb) = brdcst{1,a};
    end
    
else
    
    for i = N+1: Nc
        u = u+1;
    
        if length(brdcst{1,i})==1
            b(u,hb) = brdcst{1,i};
        elseif length(brdcst{1,i}) == 2 && i > 20
            b(u,:) = brdcst{1,i};
            hb = 2; % switch the harmonic once the length becomes 2
        end
    
    end
    
end

kk = find(b==0);
b(kk) = NaN;
fc = 25:.5:100;


%%
v = 0;
V = [];
V1 = [];


%%% Do a quick test if length(echoL{1,t}) == 20
flipE = 0; % defaulty is length(echoL{1,t}) == 20 doesn't exist


for y = N+1:Nc
    if length(echoL{1,y}) == 2*Necho
        flipE = 1;
        break;
    end
end


for t = N+1:Nc
    v = v+1;
    if length(echoL{1,t}) == Necho && t<40 % should avoid the overlap region of two harmonics, incase there is a notch, and then echoL captures 2nd harmonic here
        kkc = find(E(:,1,he));
        if ~isempty(kkc)
            TF = abs(echoL{1,t} - E(kkc(end),:,he)) < 550;
            ktf = find(TF == 0);
            if isempty(ktf)
                E(v,:,he) = echoL{1,t};
                %scatter(E(v,:,he),fc(v)*ones(10,1),'ob');
            end
            
        else
            
            E(v,:,he) = echoL{1,t};
            %scatter(E(v,:,he),fc(v)*ones(10,1),'ob');
        end

        
    elseif flipE == 1 && t > 30 && length(echoL{1,t}) == 2*Necho
        echo = echoL{1,t};
        E(v,:,1) = echo(1:2:end);
        E(v,:,2) = echo(2:2:end);
        he = 2;
        V = [V t];
        V1 = [V1 v];
        

        
        
    elseif length(echoL{1,t}) == Necho && flipE == 0 && t>67 % t>67, 67 is 53kHz, which is the end of 1st harmonic
        he = 2;
        
        kkc = find(E(:,1,he));
        if ~isempty(kkc)
            TF = abs(echoL{1,t} - E(kkc(end),:,he)) < 550;
            ktf = find(TF == 0);
            if isempty(ktf)
                E(v,:,he) = echoL{1,t};
                %scatter(E(v,:,he),fc(v)*ones(10,1),'ob');
            end
            
        else
            
            E(v,:,he) = echoL{1,t};
            %scatter(E(v,:,he),fc(v)*ones(10,1),'ob');
        end
        
    end
    
end



%%
E1 = E(:,:,1);
h1 = find(E1(:,1) == 0);
r1 = find(E1(:,1));
if length(r1)<Necho
    EF1 = R1;
    r1 = find(R1(:,1));
else
    EF1 = E1;
end



E2 = E(:,:,2);
h2 = find(E2(:,1) == 0);
r2 = find(E2(:,1));
am = intersect(h1,h2);

if length(r2)<Necho
    EF2 = R2;
    r2 = find(R2(:,1));
else
    EF2 = E2;
end

for p = 1:length(am)
    ec = echoL{1,N+am(p)};
    
    
    if am(p)<= 40 % why 40, 40+10 is 50, which corresponds to 44.5 kHz (mapped as the start Freq of 2H from spectrogram plot), 1H: 22 - 52.5kHz; 2H:47.5 - 100kHz
    % first harmonic
    [~,I] = min(abs(r1 - am(p)));
    firstHV = EF1(r1(I),:);

        for w = 1:length(ec)
            kk = find(abs(firstHV - ec(w)) < 500);
            if ~isempty(kk)
                E1(am(p),kk) = ec(w);
            end
        end
        
    elseif am(p)> 40 && am(p)<65
    
    [~,I2] = min(abs(r2 - am(p)));
    secHV = EF2(r2(I2),:);
    if ~isempty(V1)
        fHV = EF1(V1(end),:);
    else
        kj = find(EF1(:,1));
        
        fHV = EF1(kj(end),:);
    end
    
            
            for ww = 1:length(ec)
                g1 = find(abs(fHV - ec(ww))<500);
                g2 = find(abs(secHV - ec(ww))<500);
                if ~isempty(g1) && isempty(g2)
                    E1(am(p),g1) = ec(ww);
                elseif isempty(g1) && ~isempty(g2)
                    E2(am(p),g2) = ec(ww);
                elseif ~isempty(g1) && ~isempty(g2)
                    disp('error, one value was designated to two harmonics, please check\n');
                    fprintf('p value is %d',p);
                    break;
                end
            end
                    
    elseif am(p)>= 65 % 65+10 is 75, --> 57kHz, above 57kHz, only look at 2H
        
        [~,I2] = min(abs(r2 - am(p)));
        secHV = EF2(r2(I2),:);
        
        for w2 = 1:length(ec)
            kk2 = find(abs(secHV - ec(w2))<500);
            if ~isempty(kk2)
                E2(am(p),kk2) = ec(w2);
            end
        end   
            
        
    end
    
end

ee1 = find(E1 == 0);
ee2 = find(E2 == 0);
E1(ee1) = NaN;
E2(ee2) = NaN;

%%% Dechirp
D1 = E1 - b(:,1);
D2 = E2 - b(:,2);

if plotflag
figure
cc = flipud(parula(100));
des = 1:3:151;
FC = 25:.5:100;
for f = 1:Necho
    hold on, scatter(E1(des,f), FC(des), 20, cc(10*f,:),'filled','MarkerEdgeColor', cc(end,:), ...
        'MarkerFaceAlpha',.5);
    scatter(E2(des,f), FC(des), 20, cc(10*f,:),'filled','MarkerEdgeColor', cc(end,:), ...
        'MarkerFaceAlpha',.5);
end

c = jet(100);
Cc = flipud(colormap(winter(length(E1(:,1)))));
des = 1:3:151;
figure
for z = 1:Necho
    scatter(D1(des,z),fc(des), 30, 'b', 'filled', 'MarkerEdgeColor', Cc(end,:));
    hold on
    scatter(D2(des,z),fc(des), 30, 'r', 'filled', 'MarkerEdgeColor', Cc(end,:));
end

set(gca,'FontName','times','FontSize', 12);
set(gca,'XTick',0:5E3:30E3,'XTickLabel',0:10:60);
xlabel('time (ms)');
colormap('winter');
h = colorbar('eastoutside','Direction','reverse');
%ylim([25,130])
set(h, 'Ticks', [0, 1], 'TickLabels', {'high', 'low'},'FontSize', 12,'FontName', 'times');
ylabel(h, 'frequency');
end

end