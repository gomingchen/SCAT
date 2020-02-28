
function hearcellspread_empty_toomanychannels(Fc)

st_frq = Fc(1)/1E3; % in kHz
e_frq = Fc(end)/1E3; % in kHz
delta_frq = e_frq - st_frq;
figure, p = patch([0 delta_frq 0], [st_frq st_frq+delta_frq/2 e_frq], [0.39 0.58 0.93]);
p.EdgeColor = 'none';

xlabel('frequency difference (kHz)');
ylabel('frequency (kHz)');
set(gca, 'FontName', 'Times', 'FontSize', 18);
xlim([0,e_frq]);

end