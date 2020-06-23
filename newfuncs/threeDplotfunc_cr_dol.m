
function Nn_tot = threeDplotfunc_cr_dol(A,fc, q, des, step, step_thin, nbins, ipL,sep, ip_PH, lat)
% ip_ph: phantom nulls
    hold on
    minA = min(A(:,q));
    A(:,q) = A(:,q) + sep*(q-1) - minA;
    sz = 30;

    %%% Pretty scatter plot
    trans = 0.5;
    Fcs = fc; % Fc is associated with d1
    y = Fcs(des)/1E3;
    y(end) = NaN;
    Cc = flipud(colormap(winter(length(y))));
    [xs, I] = sort(A(des,q),'ascend');
    ys = y(I);
    s2 = scatter(xs,ys,sz,Cc,'filled','MarkerEdgeColor', Cc(end,:));
    s2.MarkerFaceAlpha = trans;
    th_ph = lat(q)/500 + sep*(q-1) - minA;
    plot([th_ph th_ph], [fc(1)/1E3 fc(end)/1E3], 'k','LineWidth',2);
    
    %%% Plot the histogram under the dechirped signal
    [a,edges] = histcounts(A(:,q), nbins(q));
    a1 = (edges(1)+edges(2))/2;
    a2 = edges(2) - edges(1);
    a3 = (edges(end-1) + edges(end))/2;
    bar(a1:a2:a3, a.*(fc(1)/1E3 - 5)/max(a),'FaceColor', 'c','FaceAlpha', .5, 'EdgeColor', 'black', 'LineWidth',.5);
   
    xval = max(A(:,q))+.05;
    xval_tn = xval + .2;
    mean_ip_L = ipL;

    [~, ~, C, ~, ~] = haircellspreaduniversal_upright_fig_just_triangle(fc, mean_ip_L, xval_tn);
    
    %%% Patch for the blurring effect

    blueColor2 = [0.3010 0.7450 0.9330];
    thin = 1:step_thin/1E3:length(fc);
    Fc_thin = fc(1):step_thin:fc(end);
    tf = ~isnan(ip_PH(thin));
    k_ph = find(tf);
    NL_thin = length(thin) - 1;
    Nn_ph = zeros(1, NL_thin); % count number of turned-on phantom nodes
    if ~isempty(k_ph) && length(k_ph)>5
        num = length(k_ph); % number of detection nodes turned on due to ALT
        Xs = [];Ys = []; Zs = [];
        k_ph = sort(k_ph, 'descend');
        dif_kph = -diff(k_ph);
        K = find(dif_kph>1);
        
        if isempty(K)&&~isempty(k_ph)
            K_ph = k_ph;
            
            Nk = length(K_ph);
            hold on
            scatter3((xval_tn).*ones(Nk,1), Fc_thin(K_ph)/1E3 ,zeros(Nk,1),10,'oy','filled');
            % plot the blurring effect
            
            F1 = Fc_thin(K_ph(1))/1E3;
            F2 = Fc_thin(K_ph(end))/1E3;
            X = (xval_tn).*ones(1,3);
            Y = [F1 F2+abs(F1-F2)/2 F2];
            Z = [0 (F1-F2) 0];
            p = patch(X, Y, Z, 'y', 'FaceAlpha',.3,'EdgeColor','b'); % 
            p.AlignVertexCenters = 'on';
            Xs = [Xs X]; Ys = [Ys Y]; Zs = [Zs Z];
            % for phantom nodes
            Nn_ph = Nn_ph + [(Nk-1):-1:1 zeros(1, NL_thin-(Nk-1))];
            
        elseif ~isempty(K)
            N = length(K);
            
            for a = 1:N+1
            if a == 1
                t1 = 1;
            else
                t1 = 1+K(a-1);
            end
            
            if a == N+1
                t2 = length(k_ph);
            else
                t2 = K(a);
            end
            
            K_ph = k_ph(t1:t2);
            Nk = length(K_ph);
            hold on
            scatter3((xval_tn).*ones(Nk,1), Fc_thin(K_ph)/1E3 ,zeros(Nk,1),10,'oy','filled');
            % plot the blurring effect
            if length(K_ph) > 1
            F1 = Fc_thin(K_ph(1))/1E3;
            F2 = Fc_thin(K_ph(end))/1E3;
            X = (xval_tn).*ones(1,3);
            Y = [F1 F2+abs(F1-F2)/2 F2];
            Z = [0 (F1-F2) 0];
            p = patch(X, Y, Z, 'y', 'FaceAlpha',.3, 'EdgeColor','b'); %
            p.AlignVertexCenters = 'on';
            Xs = [Xs X]; Ys = [Ys Y]; Zs = [Zs Z];
            % for phantom nodes
            Nn_ph = Nn_ph + [(Nk-1):-1:1 zeros(1, NL_thin-(Nk-1))];
            end
        
            end
            
            
        end
        
        X1 = xval_tn.*ones(1, length(Xs) + 3);
        Y1 = [fc(end)/1E3 Ys fc(1)/1E3 fc(1)/1E3+(fc(end)-fc(1))/2/1E3]; % [A ... C B]
        Z1 = [0 Zs 0 (fc(end)-fc(1))/1E3];
        Lia = ismember([fc(end)/1E3 fc(1)/1E3], Ys);
        newin = logical([1-Lia(1) ones(1,length(Ys)) 1-Lia(2) 1]); 
        p = patch(X1(newin), Y1(newin), Z1(newin), blueColor2, 'FaceAlpha',.2); % ,'LineStyle','none'
        p.EdgeColor = 'b';
        plot3(X1(1:2), Y1(1:2), Z1(1:2),'k');
        
    else
        X = xval_tn.*ones(3,1);
        Y = [fc(end)/1E3 fc(1)/1E3 fc(1)/1E3+(fc(end)-fc(1))/2/1E3];
        Z = [0 0 (fc(end)-fc(1))/1E3];
        plot3(X(1:2), Y(1:2), Z(1:2),'k');
        p = patch(X, Y, Z, blueColor2, 'FaceAlpha',.2); %,'LineStyle','none'
        p.AlignVertexCenters = 'on';
        p.EdgeColor = 'b';
        ylabel('frequency (kHz)');
    end
    
    %%% Histogram of glint spacing on xz plane. 
    
    edges = 0:step_thin/1E3:(fc(end)-fc(1))/1E3;
    [yfn,~] = histcounts(C, edges); 
    
    Nn_tot = Nn_ph + yfn;
    scale = 1/((step*1E-3)*sum(Nn_tot));
    nnp = Nn_tot.*scale;
    kn = find(nnp);
    NP = nnp(kn);
    xn = [NP;NP];
    zn = repmat([edges(kn) edges(kn(end)+1)], 2, 1);
    Xn = xval_tn - [0; xn(:); 0];
    Yn = (fc(end)/1E3).*ones(length(Xn),1);
    Zn = zn(:);

    pp = patch(Xn, Yn, Zn, 'b');
    pp.FaceAlpha = 0.3;

end