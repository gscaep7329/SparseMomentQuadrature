%
clc; clear all; close all;
randn('state',1);rand('twister',1);

rv(1) = struct('type','norm','par',[]);
rv(2) = struct('type','t','par',20);
rv(3) = struct('type','logn', 'par',[1,0.5]);
rv(4) = struct('type','exp','par',1/3);
rv(5) = struct('type','unif','par',[]);
rv(6) = struct('type','gam','par',[0.5,2]);
rv(7) = struct('type','beta','par',[1.5,0.5]);
rv(8) = struct('type','ev','par',[2,4]);
rv(9) = struct('type','weib','par',[1,5]);
rv(10) = struct('type','asin','par',[]);

dim = numel(rv);
fun = @(x) sum(x.^2,2);


%% MC simualtion
Nmc = 1e7;
x = zeros(Nmc,dim);
x(:,1) = normrnd(0, 1, Nmc,1);
x(:,2) = trnd(rv(2).par, Nmc,1);
x(:,3) = lognrnd(rv(3).par(1), rv(3).par(2), Nmc,1);
x(:,4) = exprnd(rv(4).par, Nmc,1);
x(:,5) = unifrnd(0, 1, Nmc,1);
x(:,6) = gamrnd(rv(6).par(1), rv(6).par(2), Nmc,1);
x(:,7) = betarnd(rv(7).par(1), rv(7).par(2) ,Nmc,1);
x(:,8) = ev1rnd(rv(8).par(1), rv(8).par(2), Nmc,1);
x(:,9) = wblrnd(rv(9).par(1), rv(9).par(2), Nmc,1);
x(:,10) = 2*sin(pi*rand(Nmc,1)./2).^2-1; % inverse method, see wiki of arcsin distribution


y = fun(x);
mv_mc = mean(y);
var_mc = var(y);
skew_mc = skewness(y);
kurt_mc = kurtosis(y);
fprintf(1, 'MC mean %g, var  %g, Nsim=%E\n',mv_mc,var_mc,Nmc);
fprintf(1, 'MC skew %g, kurt %g\n', skew_mc, kurt_mc);
[mv_mc,var_mc,skew_mc,kurt_mc,Nmc]


%% full tensor
nGP = 5;
M = zeros(dim,2*nGP);
W = zeros(dim, nGP);
D = zeros(dim, nGP);
for i=1:dim
    M(i,:) = raw_moment_gen(2*nGP,rv(i).type,rv(i).par);
    [~,W(i,:),D(i,:)] = aPC_Amat(nGP,M(i,:));
end

idxvec=1:nGP;
for j=2:dim
    idxvec=combvec(idxvec,1:nGP);
end
idxvec=idxvec';

zeta = zeros(size(idxvec)); weight = zeros(size(idxvec));
for i=1:size(zeta,1)
    for j=1:dim
        zeta(i,j) = D(j,idxvec(i,j));
        weight(i,j) = W(j,idxvec(i,j));
    end
end
wgt = prod(weight,2);
mv_mq = fun(zeta)'*wgt;
var_mq = (fun(zeta).^2)'*wgt-mv_mq*mv_mq;
fprintf(1,'direct moment quadrature (DMQ), mean %12.8g, var  %12.8g, Nfun=%d\n',mv_mq,var_mq,size(zeta,1));
std_mq = sqrt(var_mq);
skew_mq = (((fun(zeta)-mv_mq)/std_mq).^3)'*wgt; % https://en.wikipedia.org/wiki/Skewness
kurt_mq = (((fun(zeta)-mv_mq)/std_mq).^4)'*wgt; % https://en.wikipedia.org/wiki/Kurtosis
fprintf(1,'direct moment quadrature (DMQ), skew %12.8g, kurt %12.8g\n', skew_mq, kurt_mq);
[mv_mq,var_mq,skew_mq,kurt_mq,size(zeta,1),nGP]

%% Smolyak rule
L = 4; % Smolyak rule level,
[smol_zeta,smol_wgt]=Smolyak_rule(dim,L,rv);
mv_smq = fun(smol_zeta)'*smol_wgt(:);
var_smq = (fun(smol_zeta).^2)'*smol_wgt(:)-mv_smq*mv_smq;
fprintf(1,'sparse moment quadrature (SMQ), mean %g, var  %g, Nfun=%d\n',mv_smq,var_smq, size(smol_zeta,1));
std_smq = sqrt(var_smq);
skew_smq = (((fun(smol_zeta)-mv_smq)/std_smq).^3)'*smol_wgt(:); % https://en.wikipedia.org/wiki/Skewness
kurt_smq = (((fun(smol_zeta)-mv_smq)/std_smq).^4)'*smol_wgt(:); % https://en.wikipedia.org/wiki/Kurtosis
fprintf(1,'sparse moment quadrature (SMQ), skew %g, kurt %g\n', skew_smq, kurt_smq);
[mv_smq,var_smq,skew_smq,kurt_smq,size(smol_zeta,1),L]


exact_mean = sum(M(:,2));
jdx=combvec(1:10,1:10)';
idx=jdx(jdx(:,1)~=jdx(:,2),:);
exact_var = sum( M(idx(:,1),2).*M(idx(:,2),2) )+sum(M(:,4)) - exact_mean^2;

% THIS WORKS but above is more elegant
e2=0;
for k=1:size(jdx,1)
    kdx = jdx(k,:);
    [v,w]=unique(kdx,'stable');
    dummy = 1;
    for j=1:numel(v)
        dummy = dummy * M(v(j),2*sum(kdx==v(j)));
    end
    e2=e2+dummy;
end


% E(Y^3)=E( (sum^10_(i=1) x_i^2)^3 )
%       =E( (x_1^6+x_2^6+...+x_10^6) + 3*(x_1^4)*x_2^2+...+ )
jdx = 1:dim;
jdx = combvec(jdx,1:dim);
jdx = combvec(jdx,1:dim)';

e3=0;
for k=1:size(jdx,1)
    kdx = jdx(k,:);
    [v,w]=unique(kdx,'stable');
    dummy = 1;
    for j=1:numel(v)
        dummy = dummy * M(v(j),2*sum(kdx==v(j)));
    end
    e3=e3+dummy;
end
exact_skew = (e3-3*exact_mean*exact_var-exact_mean^3)/exact_var^1.5;

% E(Y^4)=E( (sum^10_(i=1) x_i^2)^4 );
jdx = 1:dim;
jdx = combvec(jdx,1:dim);
jdx = combvec(jdx,1:dim);
jdx = combvec(jdx,1:dim)';
e4=0;
for k=1:size(jdx,1)
    kdx = jdx(k,:);
    [v,w]=unique(kdx,'stable');
    dummy = 1;
    for j=1:numel(v)
        dummy = dummy * M(v(j),2*sum(kdx==v(j)));
    end
    e4=e4+dummy;
end
exact_kurt = (e4-4*e3*exact_mean+6*e2*exact_mean^2-3*exact_mean^4)/exact_var^2;
fprintf(1,'Exact:\n');
[exact_mean,exact_var,exact_skew,exact_kurt]