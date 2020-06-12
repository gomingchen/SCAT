% adjustment to the figure for 3D scat plot
function adjfig_3dscatCR(not)
f.Renderer = 'painters';
opengl hardware
hold off
% ax = gca;
% ax.XAxis.Visible = 'off'; % remove x-axis
set(gca, 'ZTick', [0, 5, 10:10:80], 'ZTickLabel', {'   ','200', '100', ' 50', ' 33', ' 25', ' 20', ' 17', ' 14', ' 13'});
% ax = gca;
% ax.ZDir = 'reverse';
% set(gca,'ZScale', 'log');
% zlim([2,50]);
% set(gca,'ZTick',[2 10 20 50],'ZTickLabels', {'500','100','50', '20'})


set(gca, 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:2:6*2, 'XTickLabel', 7:5:37);
xlim([-2, 14]);
view(-14,58);
grid on
set(gca, 'ZGrid', 'on', 'YGrid', 'on');
ylabel('frequency (kHz)')
xlabel('overall delay (ms)');
zlabel('glint delay ({\mu}s)');
title(sprintf('No.%d threshold', not));
