
function patchplot(edge1, edge2, h)

X = [edge1 edge1 edge2 edge2];
Y = [0 h h 0];
patch(X, Y, 'r');

end