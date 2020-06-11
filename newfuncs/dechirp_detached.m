
function [D1, D2] = dechirp_detached(brdcst, echoL, Necho)

%Necho = 10;
N = 10; % only count 25kHz and up
Nc = length(brdcst);
% for broadcast

he = 1; % harmonics flip for echoes
E = zeros(Nc-N, Necho, 2);
u = 0;
b = zeros(Nc-N,2);
hb = 1; % harmonics indicator

for i = N+1: Nc
    u = u+1;
    
    if length(brdcst{1,i})==1
        b(u,hb) = brdcst{1,i};
    elseif length(brdcst{1,i}) == 2 && i > 20
        b(u,:) = brdcst{1,i};
        hb = 2; % switch the harmonic once the length becomes 2
    end
    
end

kk = find(b==0);
b(kk) = NaN;
fc = 25:.5:100;

v = 0;

%figure
%hold on

%%% Do a quick test if length(echoL{1,t}) == 20
flip = 0; % defaulty is length(echoL{1,t}) == 20 doesn't exist
for y = N+1:Nc
    if length(echoL{1,y}) == 2*Necho
        flip = 1;
        break;
    end
end


for t = N+1:Nc
    v = v+1;
    if length(echoL{1,t}) == Necho

        E(v,:,he) = echoL{1,t};
         
    elseif flip == 1 && t > 30 && length(echoL{1,t}) == 2*Necho
        echo = echoL{1,t};
        E(v,:,1) = echo(1:Necho);
        E(v,:,2) = echo(Necho+1:end);
        he = 2;
        
    elseif flip == 0 && t>75
        he = 2;
    end
    
end




E1 = E(:,:,1);
E2 = E(:,:,2);

ee1 = find(E1 == 0);
ee2 = find(E2 == 0);
E1(ee1) = NaN;
E2(ee2) = NaN;

D1 = E1 - b(:,1);
D2 = E2 - b(:,2);

% plot figures for checking
% figure
% cc = flipud(parula(100));
% des = 1:3:151;
% FC = 25:.5:100;
% for f = 1:Necho
%     hold on, scatter(E1(des,f), FC(des), 20, 'or','filled','MarkerEdgeColor', cc(end,:), ...
%         'MarkerFaceAlpha',.5);
%     scatter(E2(des,f), FC(des), 20, 'ob','filled','MarkerEdgeColor', cc(end,:), ...
%         'MarkerFaceAlpha',.5);
% end
% 
% 
% 
% figure
% Cc = flipud(colormap(winter(length(E1(:,1)))));
% for z = 1:Necho
%     scatter(D1(des,z),fc(des), 30, 'b', 'filled', 'MarkerEdgeColor', Cc(end,:));
%     hold on
%     scatter(D2(des,z),fc(des), 30, 'r', 'filled', 'MarkerEdgeColor', Cc(end,:));
% end

end