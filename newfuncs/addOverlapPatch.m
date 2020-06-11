
%%% Add patch that indicates the overlapped area of two harmonics
% 49 kHz to 55.5 kHz
function addOverlapPatch(f1, f2, xtime)
    hold on
   
    pp = patch([-5 -5 xtime xtime], [f1 f2 f2 f1], zeros(4,1), 'g');
    pp.FaceAlpha = .3;
    pp.EdgeAlpha = 0;

end
