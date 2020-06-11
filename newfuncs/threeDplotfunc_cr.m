
function Nn_tot = threeDplotfunc_cr(A, q, des, step, step_thin, nbins, ipL,sep, ip_PH, lat)
% ip_ph: phantom nulls
    hold on
    minA = min(A(:,q));
    A(:,q) = A(:,q) + sep*(q-1) - minA;
    %scatter(A(:,q), Fc/1E3, 70, '.b');
    
    %%% pretty scatter plot
    sz = 30;
    trans = .7;
    bluecolor = [0 0.4470 0.7410];
    %des = 1:8:321;

    %%% pretty scatter plot
    trans = 0.5;
    Fcs = 20E3:.5E3:100E3; % Fc is associated with d1
    y = Fcs(des)/1E3;
    y(end) = NaN;
    Cc = flipud(colormap(winter(length(y))));
    %C = flipud([zeros(size(y))', zeros(size(y))', linspace(0,1,length(y))']);
    [xs, I] = sort(A(des,q),'ascend');
    ys = y(I);
    s2 = scatter(xs,ys,sz,Cc,'filled','MarkerEdgeColor', Cc(end,:));
    s2.MarkerFaceAlpha = trans;
    th_ph = lat(q)/500 + sep*(q-1) - minA;
    plot([th_ph th_ph], [20 100], 'k','LineWidth',2);

%     patch(xs,ys,c,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
    
    %%% Plot the histogram under the dechirped signal
    [a,edges] = histcounts(A(:,q), nbins(q));
    a1 = (edges(1)+edges(2))/2;
    a2 = edges(2) - edges(1);
    a3 = (edges(end-1) + edges(end))/2;
    bar(a1:a2:a3, a.*10/max(a),'FaceColor', 'c','FaceAlpha', .5, 'EdgeColor', 'black', 'LineWidth',.5);
   
    ylim([0, 100]);

    %[maxi, U] = max(a);
    xval = max(A(:,q))+.05;
    xval_tn = xval + .2;
     mean_ip_L = ipL;
    %[p, ~, C] = hearcellspreaduniversal_upright_fig_just_triangle(Fc,
    %mean_ip_L,1, xval_tn); % for previous sequence
    Fct = 20E3:.5E3:100E3; % leave the Fc for triangular plot uninterrupted
    NL = length(step:step:80E3);
    [p, ~, C, xc, yc] = hearcellspreaduniversal_upright_fig_just_triangle(Fct, mean_ip_L, xval_tn);
    
    %[B,~] = histc(C,unique(C));

    %%% patch for the blurring effect
%     blueColor = [0 0.4470 0.7410];
    blueColor2 = [0.3010 0.7450 0.9330];
    thin = 1:(step_thin/1E3*2):length(ip_PH); % because original step was .5kHz
    Fc_thin = 20E3:step_thin:100E3;
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
            p = patch(X, Y, Z, 'y', 'FaceAlpha',.3,'LineStyle','none');
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
            p = patch(X, Y, Z, 'y', 'FaceAlpha',.3,'LineStyle','none');
            p.AlignVertexCenters = 'on';
            Xs = [Xs X]; Ys = [Ys Y]; Zs = [Zs Z];
            % for phantom nodes
            Nn_ph = Nn_ph + [(Nk-1):-1:1 zeros(1, NL_thin-(Nk-1))];
            end
        
            end
            
            
        end
        
        
    

        % A:xval_tn, 100, 0
        % C:xval_tn, 20, 0
        % B:xval_tn, 60, 80
        X1 = xval_tn.*ones(1, length(Xs) + 3);
        Y1 = [100 Ys 20 60]; % [A ... C B]
        Z1 = [0 Zs 0 80];
        Lia = ismember([100 20], Ys);
        newin = logical([1-Lia(1) ones(1,length(Ys)) 1-Lia(2) 1]); 
        patch(X1(newin), Y1(newin), Z1(newin), blueColor2, 'FaceAlpha',.2,'LineStyle','none');
       
        
    else
        X = xval_tn.*ones(3,1);
        Y = [100 20 60];
        Z = [0 0 80];
        
        p = patch(X, Y, Z, blueColor2, 'FaceAlpha',.2,'LineStyle','none');
        p.AlignVertexCenters = 'on';
        % xlabel('frequency difference (kHz)');
        ylabel('frequency (kHz)');
    end
    %%%
    
    
    
%     if ismember(q, [5, 6])
%         [a,edges2] = histcounts(C, 8, 'Normalization', 'pdf'); 
%     else
%         [a,edges2] = histcounts(C, 'Normalization', 'pdf');
%     end

    
    %%% histogram of glint spacing on xz plane. 
    
%     if length(edges2)>2
%        edges2(1) = edges2(1)-0.5;
%        edges2(end) = edges2(end)+.5;
%     else
%         edges2(1) = edges2(1)-0.5;
%     end
    
%     for ii = 1:length(edges2)-1
%         v = a;
%         X = [xval_tn-v(ii); xval_tn; xval_tn; xval_tn-v(ii)];
%         Y = 100.*ones(4,1);
%         Z = [edges2(ii); edges2(ii); edges2(ii+1); edges2(ii+1)];
%         pr = patch(X,Y,Z,'r');
%         pr.FaceAlpha = 0.3;
%         
%     end
    
    % histogram of the distribution of total nodes
    
    % nodes from actual notches    
%     Nfn = floor(xc/(step*1E-3)); % step was in Hz
%     yfn = zeros(1,NL);
%     yfn(Nfn) = yc;
    edges = 0:step_thin/1E3:80;
    [yfn,bb] = histcounts(C, edges); 
    
    Nn_tot = Nn_ph + yfn;
    %Nn_tot = Nn_ph;
    scale = 1/((step*1E-3)*sum(Nn_tot));
    nnp = Nn_tot.*scale;
    kn = find(nnp);
    NP = nnp(kn);
    xn = [NP;NP];
    zn = repmat((0:length(kn)),2, 1)*step_thin/1E3;
    Xn = xval_tn - [0; xn(:); 0];
    Yn = 100.*ones(length(Xn),1);
    Zn = zn(:);
    pp = patch(Xn, Yn, Zn, 'b');
    pp.FaceAlpha = 0.3;
    


end