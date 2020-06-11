function [p, P, C, A,B] = haircellspreaduniversal_upright_fig_just_triangle(Fc, I, xval)
%%% This function for unequal spacing
%%% Modified by treating the dimension of Fc as a parameter
%%% step, the width of the frequency band, in kHz, for example: .5kHz or
%%% 1kHz
    N = length(Fc);
    P = zeros(length(Fc), length(Fc));
    P(:,1) = ones(length(Fc),1);
    NI = length(I);
    Nc = Fc/1E3;
%     maxF = Fc(end)/1E3;
%     minF = Fc(1)/1E3;
%     Fcm = Fc/1E3 - 1;

    %scatter3(xval.*ones(N,1), P(:,1)+ Fcm',zeros(N,1), 10, 'MarkerEdgeColor', [100/255, 149/255, 237/255]);
%     G = zeros(N,N);
%     G(:,1) = Fc/1E3;
    %step = 0.5; % write as a parameter
    
%     for g = 2:N
% 
%         hold on,
%         a = minF + step/2*(g-1);
%         b = a + (N-1-(g-1))*step;
%         G(g:end,g) = a:step:b;
%         in = g:N;
%         n = length(in);
% %         P(g:end,g) = ones(length(in),1);
% %         G(g:end,g) = P(g:end,g)+(minF-1:step:(maxF-g+.5))'+step/2*(g-1);
%         %scatter3(xval.*ones(n,1), G(g:end,g),(g-1)*step*ones(n,1),10, 'MarkerEdgeColor',[100/255, 149/255, 237/255]);
%     end
    
if length(I)>1
    

    
    %scatter(ones(NI,1), Nc(I),'*r');
    sorted_I = sort(I,'ascend');
%     diff_I = diff(sorted_I);
%     dif = diff_I(1);
%     vec = zeros(NI, NI); % y-coordinates of the spikes
%     x = zeros(1,NI); % x-coordinates of the spikes
    
    F = Nc(sorted_I);
    siz = 20;
    scatter3(xval.*ones(NI,1), F,zeros(NI,1),siz, 'or','filled');
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
    
        scatter3(xval.*ones(length(cx{2}),1), cy{2},cx{2},siz, 'or','filled');
        XX = xval.*ones(length(cx{2}),1);
        YY = cy{2};
        ZZ = cx{2};
        plot3([XX(1) XX(end)], [YY(1) YY(end)], [ZZ(1) ZZ(end)],'--','Color', [.5, .5, .5]);
        %[a,~] = histc(cx{2},unique(cx{2}));
        
        A = unique(cx{2});
        [B,~] = histc(cx{2},unique(cx{2}));
        
        
        C = cx{2};
%         bar(unique(cx{2}), a.*15/max(a),'r')
    for d = 2
        xcor = cx{d-1};
        xcor2 = cx{d};
                    
        ycor = cy{d-1};
        ycor2 = cy{d};
        for s = 1:length(xcor2)
            plot3(xval.*[1, 1], [ycor(s),ycor2(s)],[xcor(s),xcor2(s)],'-r');
            plot3(xval.*[1, 1], [ycor(s+1),ycor2(s)],[xcor(s+1), xcor2(s)], '-r');
        end
    end
    
elseif length(I) == 1
        sorted_I = sort(I,'ascend');
        F = Nc(sorted_I);
        siz = 20;
        scatter3(xval.*ones(NI,1), F,zeros(NI,1),siz, 'or','filled');
        p = [];
        C = [];
        A = [];
        B = [];
        P = [];
else
    
    p = [];
    C = [];
    A = [];
    B = [];
    P = [];
    
end

% X1 = 0;
% Y1 = Fc(1)/1E3;
% 
% X2 = 0;
% Y2 = Fc(end)/1E3;
% 
% X3 = (Fc(end)-Fc(1))/1E3;
% Y3 = (Fc(end)-Fc(1))/2/1E3 + Fc(1)/1E3;
% 
% plot3(xval.*[1, 1], [Y1, Y2], [X1, X2],'b');
% plot3(xval.*[1, 1], [Y1, Y3], [X1, X3],'b');
% plot3(xval.*[1, 1], [Y2, Y3], [X2, X3], 'b');




end