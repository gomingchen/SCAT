% adjustment to the figure for 3D scat plot
function adjfig_3dscat_dolCR(not)
f.Renderer = 'painters';
opengl hardware
hold off
set(gca, 'ZTick', [0, 10:20:110], 'ZTickLabel', {'   ', '100', ' 33', ' 20', ' 14', ' 11', '  9'});


set(gca, 'FontName', 'arial', 'FontSize', 12);
set(gca, 'XTick', 0:2:4, 'XTickLabel', '');
xlim([-1, 3]);
view(-14,58);
grid on
set(gca, 'ZGrid', 'on', 'YGrid', 'on');
xlabel('overall delay (ms)')
ylabel('frequency (kHz)');
zlabel('glint spacing ({\mu}s)');
title(sprintf('No.%d threshold', not));