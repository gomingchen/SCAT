function [firstCol, secCol] = triangularNetwork(Fc, I, cspec1, spec2)
%%% This function for unequal spacing
%%% Modified by treating the dimension of Fc as a parameter
%%% inputs: step, width of freq band, in kHz

hold on;    
if ~isempty(I)
    NI = length(I);
    Nc = Fc/1E3;
    %scatter(ones(NI,1), Nc(I),'*r');
    sorted_I = sort(I,'ascend');
    diff_I = diff(sorted_I);
    dif = diff_I(1);
    vec = zeros(NI, NI); % y-coordinates of the spikes
    x = zeros(1,NI); % x-coordinates of the spikes
    
    F = Nc(sorted_I);
    scatter(zeros(NI,1), F, ['o', cspec1], 'filled');
    firstCol = [zeros(NI,1), F'];
    xx{1} = zeros(NI,1)';
    yy{1} = F;
    
    for t = 1:NI-1
        I2 = F(t+1:end) - F(t);

%             scatter(I2, .5*(F(t)+F(t+1:end)), '*r');
        xx{t+1} = I2;
        yy{t+1} = .5*(F(t)+F(t+1:end));

    end
    NX = xx{1};
    NY = yy{1};
    
    for h = 1:NI-1
        
        nx{h} = [NX(h) xx{h+1}];
        ny{h} = [NY(h) yy{h+1}];
        if h == NI-1
            nx{h+1} = NX(h+1);
            ny{h+1} = NY(h+1);
        end
    end

    
    for q = 1:NI
        cx{q} = [];
        cy{q} = [];
    for p = 1:NI
        tx = nx{p};
        ty = ny{p};

        
            if length(tx)>=q
                cx{q} = [cx{q} tx(q)];
                cy{q} = [cy{q} ty(q)];
            end
        
    end
    end
        
        scatter(cx{2}, cy{2}, ['o', cspec1], 'filled');
        secCol = [cx{2}' cy{2}'];
%         [a,~] = histc(cx{2},unique(cx{2}));
%         xbar = unique(cx{2});
%         bar(unique(cx{2}), a.*15/max(a),'r')
    for d = 2
        xcor = cx{d-1};
        xcor2 = cx{d};
                    
        ycor = cy{d-1};
        ycor2 = cy{d};
        for s = 1:length(xcor2)
            plot([xcor(s),xcor2(s)], [ycor(s),ycor2(s)],spec2);
            plot([xcor(s+1), xcor2(s)], [ycor(s+1),ycor2(s)], spec2);
        end
    end
    
end

% xlabel('frequency difference (kHz)');
% ylabel('frequency (kHz)');
% set(gca, 'FontName', 'Times', 'FontSize', 18);
% xlim([0,80]);


end
