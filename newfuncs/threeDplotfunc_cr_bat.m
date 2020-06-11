
function threeDplotfunc_cr_bat(A, q, des, step, index,ip1, ipL,sep,Fco)
% ip_ph: phantom nulls
% Fco, Fc original
    
    minA = min(A(:,q));
    A(:,q) = A(:,q) + sep*(q-1) - minA +.2;
  
    xval = max(A(:,q))+.05;
    xval_tn = xval + .6;
    mean_ip_L = ipL;
    NL = length(step:step:80E3);
    
    [p, ~, C, xc, yc] = hearcellspreaduniversal_upright_fig_just_triangle(Fco, mean_ip_L, xval_tn);
 
    if ismember(q,index)
    Y1 = Fco(1)/1E3;
    Y2 = Fco(end)/1E3;
    Y3 = Y1 + abs(Y2-Y1)/2;
    %pp = patch([X1 X2 X3], [Y1 Y2 Y3], [Z1 Z2 Z3],'b');
    pp = patch(xval_tn*ones(3,1), [25, 100, Y3], [0 0 75],'b');
    pp.FaceAlpha = .1;
    pp.EdgeColor = 'b';
    pp.EdgeAlpha = .5;
    end
    
    % histogram of the distribution of total nodes

 if ~isempty(C) 
     c = unique(C);
        if length(C) == 1
            edges2 = [C-1, C+1];
            a = 1;
        else
            [a,edges2] = histcounts(C,length(c));
        end
        
    if length(ip1)>=2
        
        pacolor = 'r';
    else
        pacolor = 'b';
    end


    [~,Im] = max(a);
    if exist('edges2','var') && ismember(q, index)
        str = sprintf('%d \\mus', round(1000/((edges2(Im)+edges2(Im+1))/2)));
        
    end
    
    if length(edges2)>2
       edges2(1) = edges2(1)-0.5;
       edges2(end) = edges2(end)+.5;
%         V = edges2(1);
%         delta = edges2(2) - edges2(1);
%         n = 3;
%         if rem(length(edges2), 2) == 0 % even
%             b = V + (length(edges2)/2 - 1)*(1-n)*delta;
%         else
%             b = V + (length(edges2)-1)/2*(1-n)*delta;
%         end
%         newE = b:n*delta:b+(length(edges2)-1)*n*delta;
%         edges2 = newE;
%     else
%         edges2(1) = edges2(1)-0.5;
    end
    

    for ii = 1:length(edges2)-1
        v = 0.7*a/max(a);
        X = [xval_tn-v(ii); xval_tn; xval_tn; xval_tn-v(ii)];
        Y = 100.*ones(4,1);
        Z = [edges2(ii); edges2(ii); edges2(ii+1); edges2(ii+1)];
        patch(X,Y,Z,pacolor);
    end
    
    if exist('str','var')
        text(xval_tn-1, 100, Z(1)+10, str, 'FontSize', 12, 'FontName', 'times');
    end
    
 end
    



end
