

%addpath('figmats');
%load('eighttargets.mat');

[row, col] = find(vec_all);
browncolor = [56 28 0]/255;
greycolor = [152, 164, 174]/255;
darkgrey = [1 1 1 ]*.5;
cc = colormap(hsv(tars.Ns));
cb = colormap(winter(tars.Ns));
gs = round(gs);
cc(3,:) = cb(3,:);
figure
hold on
for i = 1:tars.Ns
    ind = i;
    coord_ss = [tars.r(ind)*cosd(tars.theta(ind)+90), tars.r(ind)*sind(tars.theta(ind)+90)];
    scatter(coord_ss(1), coord_ss(2), 30, 'o', 'MarkerEdgeColor', cc(i,:), 'MarkerFaceColor', cc(i,:));
    
    text(coord_ss(1), coord_ss(2)+.6, sprintf('%d',i), 'FontName','times', 'FontSize', 12);
    
    %%% scenario 1
%     if ismember(i, [1 2 4 6])
%         text(coord_ss(1)-2.2, coord_ss(2), [sprintf('(%d ', gs(i)) '{\mu}s)'],'Color',darkgrey, 'FontName', 'Times', 'FontSize', 10);
%     elseif ismember(i,[3 5])
%         text(coord_ss(1)+.3, coord_ss(2), [sprintf('(%d ', gs(i)) '{\mu}s)'], 'Color', darkgrey, 'FontName', 'Times', 'FontSize', 10);
%     elseif ismember(i,8)
%         text(coord_ss(1)-2.2, coord_ss(2), [sprintf('(%d ', gs(i-1)) '{\mu}s)'],'Color',darkgrey, 'FontName', 'Times', 'FontSize', 10);
%     end
    
    %%% scenario 2
%     if ismember(i, [1 2 4])
%         text(coord_ss(1)-2.2, coord_ss(2), [sprintf('(%d ', gs(i)) '{\mu}s)'],'Color',darkgrey, 'FontName', 'Times', 'FontSize', 10);
%     elseif ismember(i,3)
%         text(coord_ss(1)+.3, coord_ss(2), [sprintf('(%d ', gs(i)) '{\mu}s)'], 'Color', darkgrey, 'FontName', 'Times', 'FontSize', 10);
%     end
end
tv = diff(tarvec);
kt = [find(tv), length(tarvec)-1];
for e = 1:max(row)
    
    xm = (coordear_all(e,1) + coordear_all(e,3))/2;
    ym = (coordear_all(e,2) + coordear_all(e,4))/2;
    scatter([coordear_all(e,1), coordear_all(e,3)], [coordear_all(e,2), coordear_all(e,4)], 10, 'o', 'MarkerFaceColor',browncolor, 'MarkerEdgeColor', browncolor);
    plot([coordear_all(e,1), coordear_all(e,3)], [coordear_all(e,2), coordear_all(e,4)],'-b');
    plot([xm, xm+vec_all(e,1)/2], [ym, ym+vec_all(e,2)/2], 'Color', cc(tarvec(e),:));
    if ismember(e,kt)
        p2 = 30*[coordear_all(e,1)-coordear_all(e,3), coordear_all(e,2) - coordear_all(e,4)];
        %dp = [vec_all(e,1)/2+vec_all(e,2)/2, vec_all(e,2)/2-vec_all(e,1)/2];
%         h2 = quiver(xm+.2*p2(1), ym+.2*p2(2), p2(1),p2(2),'LineWidth',2,'Color',  cc(tarvec(e),:));
%         set(h2,'AutoScale','on', 'AutoScaleFactor',2)
        
        ind = tarvec(e);
        coord_ss = [tars.r(ind)*cosd(tars.theta(ind)+90), tars.r(ind)*sind(tars.theta(ind)+90)];
        plot([xm, coord_ss(1)], [ym, coord_ss(2)],':', 'Color', cc(ind,:), 'LineWidth',2);
    end
        
end
%scatter(coord_ss(1), coord_ss(2), 'or', 'filled');
axis equal
set(gca,'FontName','times', 'FontSize', 12);


xlabel('{\it x} (m)');
ylabel('{\it y} (m)');


%% inset
% figure
% Fc = 20:1:100;
% for j = 1:4
%     subplot(4,1,j);
%     simple_triNet(Fc, mean_ip{j}, 10, cc(6,:));
% end
