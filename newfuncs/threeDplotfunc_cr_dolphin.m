
function threeDplotfunc_cr_dolphin(A, q, des, nbins, ipL,sep,Fco)
% ip_ph: phantom nulls
% Fco, Fc original
    hold on
    minA = min(A(:,q));
    A(:,q) = A(:,q) + sep*(q-1) - minA;
    sz = 30;


    %%% Pretty scatter plot
    trans = 0.5;
    y = Fco(des)/1E3;
    y(end) = NaN;
    Cc = flipud(colormap(winter(length(y))));
    [xs, I] = sort(A(des,q),'ascend');
    ys = y(I);
    s2 = scatter(xs,ys,sz,Cc,'filled','MarkerEdgeColor', Cc(end,:));
    s2.MarkerFaceAlpha = trans;
    
    %%% Plot the histogram under the dechirped signal
    [a,edges] = histcounts(A(:,q), nbins(q));

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
        yy = [0 a(jj).*20/max(a) a(jj)*20/max(a) 0];
        patch(xx,yy,'c');
        
    end

    xval = max(A(:,q))+.05;
    xval_tn = xval + .2;
    mean_ip_L = ipL;

    
    [~, ~, C, ~, ~] = haircellspreaduniversal_upright_fig_just_triangle(Fco, mean_ip_L, xval_tn);
 


    Y1 = Fco(1)/1E3;
    Y2 = Fco(end)/1E3;
    Y3 = Y1 + abs(Y2-Y1)/2;
    pp = patch(xval_tn*ones(3,1), [Y1, Y2, Y3], [0 0 Y2-Y1],'b');
    pp.FaceAlpha = .1;
   
    
    % Histogram of the distribution of total nodes

    if ~isempty(C)
    [a,edges2] = histcounts(C);

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
        Y = Y2.*ones(4,1);
        Z = [edges2(ii); edges2(ii); edges2(ii+1); edges2(ii+1)];
        patch(X,Y,Z,'r');
    end
    end
    
    if exist('edges2')
    str = sprintf('%d \\mus', round(1000/edges2(Im)));
    text(xval_tn-1, Y2, Z(1)+10, str, 'FontSize', 12, 'FontName', 'arial');
    end


end