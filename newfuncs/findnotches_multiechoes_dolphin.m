
function ipL = findnotches_multiechoes_dolphin(Necho, D)
% find the notches of multiple echoes from the dechirped representations

ipL = {};
Nthresh = size(D(1).data, 2);
for p = 1:Necho
    T = D(p).data;
    [row, col] = find(isnan(T));
    [C, ia, ~] = unique(col);
    
    % start with the column (threshold) that has more than one cluster of
    % Nans
    mn = [];
    v = 0;
    while isempty(mn)
        v = v + 1;
        if (v+1)>length(ia)
            break;
        else
        row1 = row(ia(v):(ia(v+1)-1));
        mn = find(diff(row1)>1);
        end
    end
    
    if isempty(mn)
        ipL{p} = [];
        break;
    else
        
    colstart = col(ia(v));
    end
    
  
    %%%
    
    IP = [];
    ip = NaN;
    u = 0;
    
    if colstart+2<10
        G = 3;
    else
        G = 10-colstart+1;
    end
    
    for y = 1:G
    ip = findnotches(T, colstart+u);
    u = u + 1;
    if ~isnan(ip)
        IP = [IP ip'];
    end
    end
    
    
    %%%
    
 
    
    P = unique(IP);
    sort_ip = sort(P,'ascend');
    pu = 1;
    v = 5;
    difP = diff(sort_ip);
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
    warning('off','all');
    M = mode(difP);
    if (con_edges(4) - con_edges(1))<7 || min(difP)>.8*M
        v = 0;
    end
    end
    
    % eliminate ajacent notches
    kk = find(difP == 1);
    if ~isempty(kk)
        sort_ip(kk) = [];
        ipL{p} = sort_ip;
    else
        ipL{p} = sort_ip;
    end
end


end
