% adjustment to the figure for 3D scat plot
function adjfig_3dscat_dolphin(not)
f.Renderer = 'painters';
opengl hardware
hold off
% ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
%set(gca, 'ZTick', [0, 5, 10:10:80], 'ZTickLabel', {'   ','200', '100', ' 50', ' 33', ' 25', ' 20', ' 17', ' 14', ' 13'});
% ax = gca;
% ax.ZDir = 'reverse';
% set(gca,'ZScale', 'log');
% zlim([2,50]);
%set(gca,'ZTick',[5 10 20 50 100 125],'ZTickLabels', {'200','100','50', '20', '10', '8'})
set(gca, 'ZTick', []);

set(gca, 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:2:20, 'XTickLabel', 6:5:(6+5*10));
xlim([-2, 20]);
view(-14,58);
grid on
set(gca, 'ZGrid', 'on', 'YGrid', 'on');
ylabel('frequency (kHz)')
xlabel('overall delay (ms)');
zlabel('gline delay ({\mu}s)');
title(sprintf('No.%d threshold', not));