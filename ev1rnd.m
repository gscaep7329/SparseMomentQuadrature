function [r] = ev1rnd(mu,sigma,varargin)
% generate random numbers from extreme value type I distribution
if nargin < 2
    error(message('ev1rnd:TooFewInputs'));
end

sizeOut = varargin{1};

% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;
r = mu - sigma*log(-log(rand(1,sizeOut)));