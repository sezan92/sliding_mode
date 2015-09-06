function out = implicitTwisting(in)
% Explicit twisting
sigma = in(1);
sigmaDot = in(2);
K = in(3);
beta = in(4);
h = in(5);
oldLambda(1) = in(6);
oldLambda(2) = in(7);

% this could be optimized and computed only once
mat = zeros(2, 2);

mat(2,:) = [h h*beta];
mat(2,:) = K*mat(2,:);
mat(1,:) = h/2*mat(2,:);

q = [sigma + h*sigmaDot; sigmaDot];

lambda = enumTwistingFastStart(mat, q, oldLambda, h);

%
out = lambda;