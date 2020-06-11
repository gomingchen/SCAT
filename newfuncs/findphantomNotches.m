
function ip =  findphantomNotches(NE, T, ipL, Nch, Nthre)


for j = 1:Nthre
    sn = NaN(Nch, NE);
    c = colormap('winter');
    bluecolor = [0 0.4470 0.7410];
    fc = 20:.5:100;
    den = 1:3:161;
    lat_th = zeros(NE,1);
    figure
    hold on
    Nbin = [3,3,3, 4, 5,5,5];
    for i = 1:NE
        R = T(i).data;
        mnd = min(R(den,j));
        s = scatter(R(den,j)-mnd+1E3*i, fc(den), 20, 'MarkerEdgeColor', [0 .5 .5], 'MarkerFaceColor', [0 .7 .7]);
        s.MarkerFaceAlpha = 1;
        %xlim([0,22E3]);
        ind = ipL{i};
        ki = find(ind<80);
        index = ind(ki);
        Max1 = max(R(index,j));
        % define the maximum value in histogram as the threshold to turn on
        % the detection layers
        [N,edges] = histcounts(R(:,j),10);
        [N2,edges2] = histcounts(R(:,j),Nbin(i)); % for plot purposes
        [~,I] = max(N2);
        %Max = floor((edges(I) + edges(I+1))/2);
        Max = edges2(I+1);
        
        
        lat_th(i) = Max;
        %%%%%
        plot([Max Max]-mnd+1E3*i, [20, 100], 'k', 'LineWidth',1);
        set(gca, 'FontName', 'times', 'FontSize', 12);
        %set(gca, 'XTick', 0:2500:22500, 'XTickLabel',0:5:41);
        set(gca,'XTick', 1000:1000:8000,'XTickLabel', 7:5:(7+5*6));
        ylabel('frequency (kHz)');
        xlabel('time (ms)');
        % plot the distribution to check the latency threshold
        bar(edges2(1:end-1)+(edges2(2)-edges2(1))/2-mnd+1E3*i, N2/max(N2)*15, 'FaceColor','b','EdgeColor','b');
        scatter(zeros(1,length(R(den,j))), fc(den)', 20, 'MarkerEdgeColor', [0 .5 .5], 'MarkerFaceColor', [0 .7 .7]);
        kk = find(R(:,j)>Max);
        kk_nan = find(isnan(R(:,j)));
        if ~isempty(kk_nan)
            kk = sort([kk; kk_nan],'ascend');
        end
        
        for t = 1:length(index)
            
            gk = find(abs(kk - index(t))<5);
            kk(gk) = [];
            
            
        end
        
        sn(kk,i) = 1;
        
      
    end
    ip(j).ft = sn;
    ip(j).lat = lat_th;
end

end