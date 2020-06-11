
function [E1, E2] = create2Href(brdcst, echoL, Necho, echolength)

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


v = 0;
V = [];
V1 = [];
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
        kkc = find(E(:,1,he));
        if ~isempty(kkc)
            TF = abs(echoL{1,t} - E(kkc(end),:,he)) < echolength;
            ktf = find(TF == 0);
            if isempty(ktf)
                E(v,:,he) = echoL{1,t};
            end
            
        else
            E(v,:,he) = echoL{1,t};
            
        end
        
        
    elseif flip == 1 && t > 30 && length(echoL{1,t}) == 2*Necho
        echo = echoL{1,t};
        E(v,:,1) = echo(1:2:end);
        E(v,:,2) = echo(2:2:end);
        he = 2;
        V = [V t];
        V1 = [V1 v];
        
    elseif flip == 0 && t>75
        he = 2;
    end
    
end



%%
E1 = E(:,:,1);
E2 = E(:,:,2);



end