

function [W, d1] = cell2doublematrix_glint(Mat,Necho, Nthre, Nfc)
% inputs: Necho, no. of echoes following the broadcast
%         Nthre, no. of thresholds
%         Nfc, no. of frequency channels
% outputs: d1, matrix contains dechirped representations of every echo



d1 = NaN(Nfc,Necho,Nthre);
b = Mat(1,100);
base = b{1,1};
for i = 1:Nthre
    cl = Mat(i,:);
    for j = 1:length(cl)
        
        G = cl(j);
        T = G{1,1};

            for t = 1:length(T)
                kk = find(abs(base - T(t))<30);
                d1(j,kk,i) = T(t);
            end
      
        
        
    end
   
end

data = zeros(Nfc, Nthre);
for p = 1:Necho
    for r = 1:size(d1,3)
       data(:,r) = d1(:,p,r); 
    end
    W(p).data = data;
end


end