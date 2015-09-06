function out = explicitTwisting(in)
% Explicit twisting
sigma = in(1);
sigmaDot = in(2);
K = in(3);
beta = in(4);

%
u=-K*(sign(sigma)+beta*sign(sigmaDot));
%
out = u;