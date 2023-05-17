function [C,cnt] = combnk(d,k)
% Generates all d-tuples whose sum of elements is equal to k
% Input:
%   d: dimension of tuples
%   k: sum of tuple elements
% Output:
%   C: matrix of d-tuples whose sum of elements is equal to k

C = [];
if d == 1
    C = k;
elseif k == 1
    C = ones(1,d);
else
    for i = 1:(k-d+1)
        T = combnk(d-1,k-i);
        T(:,end+1) = i;
        C = [C; T];
    end
end
cnt = size(C,1);
end