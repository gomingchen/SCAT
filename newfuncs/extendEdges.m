
function [st, A] = extendEdges(edges, times, midpoint)
    
N = length(edges);
a = edges(2) - edges(1);
A = times*a;


if rem(N, 2) % ==1, odd
    %midpoint = edges(round(N/2));
    st = midpoint - floor((N-1)/2)*A;
else
    %midpoint = (edges(N/2) + edges(N/2+1))/2;
    st = midpoint - (N-1)/2*A;
end

end