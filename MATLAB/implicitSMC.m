function out = implicitSMC(in)
% computes the solution of the scalar generalized equation
% 0 \in s_k - hGv + N_[-1,1](v)
% where N_[-1,1](v) is the normal cone of the hypercube of length G at v
% Since the equation is scalar, the solution is a simple projection
%
% author: Olivier Huber <olivier.huber@inria.fr>

G = in(1);
h = in(2);
s_k = in(3);

toProj = -s_k/(G*h);

if toProj > 1
    out = G;
elseif  toProj < -1
    out = -G;
else
    out = G*toProj;
end
    