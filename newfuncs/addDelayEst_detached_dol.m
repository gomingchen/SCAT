
function t = addDelayEst_detached_dol(H1, H2,sep)
% Delay estimation of the echoes from bat pulses is decided by both of the
% harmonics, but the overlapped region of the two harmonics will not be
% counted

% 
% minA = min(H1);
% H1 = H1 + sep*(q-1) - minA;
% H2 = H2 + sep*(q-1) - minA;


H1 = H1 + sep;
H2 = H2 + sep;

hold on,
H = NaN(size(H1));
H(1:49) = H1(1:49);
H(62:end) = H2(62:end);

[a,edges] = histcounts(H, 5);

[~,Im] = max(a);
if exist('edges','var')
   t = (edges(Im)+edges(Im+1))/2;
end

V = edges(1);
delta = edges(2) - edges(1);
n = 3;
if rem(length(edges), 2) == 0 % even
        b = V + (length(edges)/2 - 1)*(1-n)*delta;
else
        b = V + (length(edges)-1)/2*(1-n)*delta;
end
    
newE = b:n*delta:b+(length(edges)-1)*n*delta;

for jj = 1:length(newE)-1
    xx = [newE(jj) newE(jj) newE(jj+1) newE(jj+1)];
    yy = [0 a(jj).*20/max(a) a(jj)*20/max(a) 0];
    patch(xx,yy,'c');
        
end


end