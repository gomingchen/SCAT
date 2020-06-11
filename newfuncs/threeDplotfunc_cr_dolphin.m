
function threeDplotfunc_cr_dolphin(A, q, des, step, step_thin, nbins, ipL,sep,Fco)
% ip_ph: phantom nulls
% Fco, Fc original
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
    %Fcs = 20E3:.5E3:100E3; % Fc is associated with d1
    y = Fco(des)/1E3;
    y(end) = NaN;
    Cc = flipud(colormap(winter(length(y))));
    %C = flipud([zeros(size(y))', zeros(size(y))', linspace(0,1,length(y))']);
    [xs, I] = sort(A(des,q),'ascend');
    ys = y(I);
    s2 = scatter(xs,ys,sz,Cc,'filled','MarkerEdgeColor', Cc(end,:));
    s2.MarkerFaceAlpha = trans;
%     th_ph = lat(q)/500 + sep*(q-1) - minA;
%     plot([th_ph th_ph], [20 100], 'k','LineWidth',2);

%     patch(xs,ys,c,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
    
    %%% Plot the histogram under the dechirped signal
    [a,edges] = histcounts(A(:,q), nbins(q));
%     a1 = (edges(1)+edges(2))/2;
%     a2 = edges(2) - edges(1);
%     a3 = (edges(end-1) + edges(end))/2;
%     bar(a1:a2:a3, a.*10/max(a),'FaceColor', 'c','FaceAlpha', .5, 'EdgeColor', 'black', 'LineWidth',.5);
    V = edges(1);
    delta = edges(2) - edges(1);
    n = 5;
    if rem(length(edges), 2) == 0 % even
        b = V + (length(edges)/2 - 1)*(1-n)*delta;
    else
        b = V + (length(edges)-1)/2*(1-n)*delta;
    end
    
    newE = b:n*delta:b+(length(edges)-1)*n*delta;

    for jj = 1:length(newE)-1
        xx = [newE(jj) newE(jj) newE(jj+1) newE(jj+1)];
        yy = [0 a(jj).*10/max(a) a(jj)*10/max(a) 0] - 20;
        patch(xx,yy,'c');
        
    end
   
    ylim([-20, 200]);

    %[maxi, U] = max(a);
    xval = max(A(:,q))+.05;
    xval_tn = xval + .2;
     mean_ip_L = ipL;
    %[p, ~, C] = hearcellspreaduniversal_upright_fig_just_triangle(Fc,
    %mean_ip_L,1, xval_tn); % for previous sequence
    %Fct = 20E3:.5E3:100E3; % leave the Fc for triangular plot uninterrupted
    NL = length(step:step:80E3);
    
    [p, ~, C, xc, yc] = hearcellspreaduniversal_upright_fig_just_triangle(Fco, mean_ip_L, xval_tn);
 

    %[B,~] = histc(C,unique(C));


    Y1 = Fco(1)/1E3;
    Y2 = Fco(end)/1E3;
    Y3 = Y1 + abs(Y2-Y1)/2;
    %pp = patch([X1 X2 X3], [Y1 Y2 Y3], [Z1 Z2 Z3],'b');
    pp = patch(xval_tn*ones(3,1), [0, 200, Y3], [0 0 200],'b');
    pp.FaceAlpha = .1;
   
    
    % histogram of the distribution of total nodes

    if ~isempty(C)
    [a,edges2] = histcounts(C);

%     xval = xval+0.4;
    [~,Im] = max(a);
    
    if length(edges2)>2
       edges2(1) = edges2(1)-0.5;
       edges2(end) = edges2(end)+.5;
    else
        edges2(1) = edges2(1)-0.5;
    end
    
    for ii = 1:length(edges2)-1
        v = 0.7*a/max(a);
        X = [xval_tn-v(ii); xval_tn; xval_tn; xval_tn-v(ii)];
        Y = 200.*ones(4,1);
        Z = [edges2(ii); edges2(ii); edges2(ii+1); edges2(ii+1)];
        patch(X,Y,Z,'r');
    end
    end
    
    if exist('edges2')
    str = sprintf('%d \\mus', round(1000/edges2(Im)));
    text(xval_tn-1, 200, Z(1)+10, str, 'FontSize', 12, 'FontName', 'times');
    end


end