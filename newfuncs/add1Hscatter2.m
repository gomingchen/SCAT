
function add1Hscatter2(A, des,sep,Fco, Cc, trans)
% add scatter plot of the 1st harmonics
% inputs, A, the dechirped representations of each frequency channel
%         q, which threshold
%         des, density, usually every three points to avoid overcrowded
%         plots
%         sep, separation between echoes
%         Fco, original frequency channels, Fc in wav parameters
%         trans, transparency of the plot
   
    
    A = A + sep;
   

    sz = 30;

    %trans = 0.8;
    y = Fco(des)/1E3;
    y(end) = NaN;
    %Cc = flipud(colormap(winter(3*length(y))));
    [xs, I] = sort(A(des),'ascend');
    ys = y(I); 
    hold on
    s1 = scatter(xs,ys,sz,Cc(1:length(y),:),'filled','MarkerEdgeColor',Cc(end,:)); % Cc(1,:)
    s1.MarkerFaceAlpha = trans;
    

end
