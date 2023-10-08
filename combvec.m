function out = combvec(A,varargin)
% Generate all combinations of input arrays
npar = numel(varargin);
if npar == 0
    out = A;
elseif npar > 1 % recursive
    out = combvec(combvec(A,varargin{1}),varargin{2:end});
else
    B = varargin{1};
    out = [repmat(A,1, numel(B));
       reshape(repmat(B,size(A,2),1), 1, size(A,2)*numel(B))];
end
   
