
function ipL = findnotches_multiechoes_bat1H(Necho, index, D, ran)
% find the notches of multiple echoes from the dechirped representations
% inputs, ran, the range of 1st harmonic derived from the broadcast
ipL = {};

for q = 1:length(index)
    p = index(q);
    G = D(p).data;
    T = G(ran(1):ran(2),:);

%%%     First, check 1st column (threshold) to see if there are notches at
%%%     the beginning or the end of 1st harmonic, and also in the middle
    IP = [];
  for colstart = 1:4
    gp = isnan(T(:,colstart));
    kgp = find(gp);
    dif_kgp = diff(kgp); % see if it's one cluster or a few
    ng = find(dif_kgp>1);
    
    if ~isempty(kgp)
        if isempty(ng) % multiple clusters of NaNs
            tem = round(mean(kgp));
            IP = [IP tem];
        elseif length(ng) == 1
                a1 = 1;
                a2 = ng;
                in1 = round(mean(kgp(a1:a2)));
                in2 = round(mean(kgp(a2+1:end)));
                IP = [IP in1 in2];
                
        elseif length(ng)>1
            for a = 1:length(ng)+1
                if a == 1
                    a1 = 1;
                    a2 = ng(a);
                elseif a == length(ng)+1
                    a1 = ng(a-1)+1;
                    a2 = length(kgp);
                else
                    a1 = ng(a-1)+1;
                    a2 = ng(a);
                end

                in = round(mean(kgp(a1:a2)));
                IP = [IP in];

            end
        end

    
    end
    
  end
    
% check for repeated notches over thresholds
    P = unique(IP);
    sort_ip = sort(P,'ascend');
    pu = 1;
    v = 5;
    difP = diff(sort_ip);
    
    kk_ahead = find(difP == 1);
    if ~isempty(kk_ahead)
        sort_ip(kk_ahead) = [];
        difP = diff(sort_ip);
    end
    
    [con, con_edges] = histcounts(difP, 3);
    [~,I] = max(con);
    maj = find(difP> con_edges(I));
    M = mode(difP(maj));
    if min(difP)<.5*M
        v = 1;
    else
        v = 0;
    end
    
    while v

        pu = find(difP<con_edges(2));
    if length(pu) == 2 && pu(2) - pu(1) == 1
        sort_ip(pu(2)) = [];
    else
        sort_ip(pu) = [];
    end
    difP = diff(sort_ip);
    [con, con_edges] = histcounts(difP, 3);
    M = mode(difP);
    if (con_edges(4) - con_edges(1))<7 || min(difP)>.8*M
        v = 0;
    end
    end

    % eliminate ajacent notches
    kk = find(difP <5);
    if ~isempty(kk)
        sort_ip(kk+1) = [];
        ipL{p} = sort_ip;
    else
        ipL{p} = sort_ip;
    end
end

leftover = setdiff(1:Necho,index);

for h = 1:length(leftover)
    c = leftover(h);
    ipL{c} = [];
end


end
