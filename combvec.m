function out = my_combvec(A,B)
% Generate all combinations of input arrays
out = [repmat(A,1, numel(B));
       reshape(repmat(B,size(A,2),1), 1, size(A,2)*numel(B))];
   
