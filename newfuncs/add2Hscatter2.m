
function add2Hscatter2(A,des,sep,Fco, Cc, trans)
% add scatter plot of the 2nd harmonics
% inputs, A, the dechirped representations of each frequency channel
%         q, which threshold
%         des, density, usually every three points to avoid overcrowded
%         plots
%         sep, separation between echoes
%         Fco, original frequency channels, Fc in wav parameters
%         trans, transparency of the plot
   
    
    A = A + sep + .2;
    sz = 30;
    y = Fco(des)/1E3;
    y(end) = NaN;
    %Cc = flipud(colormap(winter(3*length(y))));
    [xs, I] = sort(A(des),'ascend');
    ys = y(I); 
    hold on
    s2 = scatter(xs,ys,sz,Cc(2*length(y)+1:3*length(y),:),'filled','MarkerEdgeColor',Cc(round(2.5*length(y)),:)); % Cc(1,:)
    s2.MarkerFaceAlpha = trans;
    

end
