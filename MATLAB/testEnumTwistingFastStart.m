

h = 0.1;
beta = .5;
gain = 2;

nIter = 1000;

Sigma = zeros(nIter, 2);
lambda = zeros(nIter, 2);


Astar = [1, h; 0 1];

Bstar(2,:) = [h h*beta];
Bstar(2,:) = gain*Bstar(2,:);
Bstar(1,:) = h/2*Bstar(2,:);

Sigma(1, :) = [1.1, 1.1];
lambda(1, :) = [1, 1];
inputVector = [Sigma(1, :), gain, beta, h, lambda(1, :)];

for i=1:nIter
    lambda(i, :) = implicitTwisting(inputVector);
    Sigma(i+1, :) = Astar*Sigma(i, :)' + Bstar*lambda(i, :)';
    inputVector(1:2) = Sigma(i+1, :);
    inputVector(6:7) = lambda(i, :);
end

figure
plot(Sigma(:,1), Sigma(:,2))