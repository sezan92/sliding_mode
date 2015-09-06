function [lambda, nbIter] = enumTwistingFastStart(mat, q, oldLambda, h)
% Find lambda
% Sigma_{k+1} = q + mat*lambda 
% lambda in -Sgn(Sigma_{k+1})

nbIter = 0;
% we suppose here that we have a symmetric control set
Kvert = [1, 1];
posIndxOld = [];
nposIndxOld = 0;
% first try the oldvalue
lambda = oldLambda;
if abs(lambda(1)) ~= Kvert(1)
    lambda(1) = Kvert(1);
end
if abs(lambda(2)) ~= Kvert(2)
    lambda(2) = Kvert(2);
end

extremaControl = mat(2,1) + mat(2,2);
indxZero = 0;

nbIterMax = 9; %2^3 +1

alreadyDone = zeros(9, 1);

alreadyDone(lambda(1) + 2*lambda(2)+6) = 1;

while nbIter < nbIterMax 
    % see if we can reach the origin
    if (lambda(1) == 0) && (lambda(2) == 0)
        valU = -q;
        if (abs(valU(1)/valU(2) - h/2) < eps) && (abs(valU(2)) < extremaControl + eps)
            lambda(1) = valU(2)/extremaControl;
            lambda(2) = valU(2)/extremaControl;
            keyboard
            return
        else
            lambda = oldLambda2;
            % turn right
            posIndx = abs(lambda(1) + lambda(2))/2+1;
            if (posIndx ~=1) && (posIndx ~= 2)
                keyboard;
            end
            nposIndx = 1;
            posIndxOld = [];
            nposIndxOld = 0;
        end
    else
        if indxZero > 0
            % we postulated that lambda(1) is in the multivalued part
            if indxZero == 1
                valU1 = -(q(1) + mat(1, 2)*lambda(2));
                if abs(valU1) < mat(1, 1)
                    lambda(1) = valU1/mat(1, 1);
                    Sigma = q + mat*lambda';
                    almostZero = abs(Sigma) < eps;
                    Sigma(almostZero) = 0;
                else
%                    keyboard;
                end
                if sign(Sigma(2)) == -sign(lambda(2))
                    return;
                else
%                    keyboard;
                    lambda = oldLambda2;
                    posIndx = abs(lambda(1) + lambda(2))/2+1;
                    if (posIndx ~=1) && (posIndx ~= 2)
                        keyboard;
                    end
                    nposIndx = 1;
                    posIndxOld = [];
                    nposIndxOld = 0;
                end
            % we postulated that lambda(2) is in the multivalued part
            else
                valU2 = -(q(2) + mat(2, 1)*lambda(1));
                if abs(valU2) < mat(2, 2)
                    lambda(2) = valU2/mat(2, 2);
                    Sigma = q + mat*lambda';
                    almostZero = abs(Sigma) < eps;
                    Sigma(almostZero) = 0;
                else
%                    keyboard;
                end
                if sign(Sigma(1)) == -sign(lambda(1))
                    return;
                else
%                    keyboard;
                    lambda = oldLambda2;
                    posIndx = abs(lambda(1) + lambda(2))/2+1;
                    if (posIndx ~=1) && (posIndx ~= 2)
                        keyboard;
                    end
                    nposIndx = 1;
                    posIndxOld = [];
                    nposIndxOld = 0;
                end
            end
        else % lambda is one of the vertex of K
            Sigma = q + mat*lambda';
            almostZero = abs(Sigma) < eps;
            Sigma(almostZero) = 0;
            sLambda = sign(lambda)';
            sProd = sLambda.*sign(Sigma);
            
            posIndx = [];
            if sProd(1)>0
                posIndx = [posIndx, 1];
            end
            if sProd(2)>0
                posIndx = [posIndx, 2];
            end
            nposIndx = numel(posIndx);
            % if this is true, we are done
            if nposIndx == 0
                return;
            end
        end
    end
    
    % prepare next iteration
    indxZero = 0;
    oldLambda2 = lambda;
    hasZero = 0;
    for i=1:nposIndx
        indx = posIndx(i);
        if (nposIndxOld == nposIndx) && ((( nposIndxOld >= 1 ) && (indx == posIndxOld(1))) || (( nposIndxOld >= 2 ) && (indx == posIndxOld(2))))
            lambda(indx) = 0;
            indxZero = indx;
            hasZero = 1;
        elseif hasZero == 0
            lambda(indx) = lambda(indx) - 2*sign(lambda(indx));
            % normalize
            % XXX does not work when K is not the unit hypercube
            lambda(indx) = lambda(indx)/abs(lambda(indx));
        end
    end
    nbIter = nbIter + 1;
    posIndxOld = posIndx;
    nposIndxOld = nposIndx;
     
    if hasZero == 0 && (lambda(2) ~= 0)
        indxConfig = lambda(1) + 2*lambda(2) + 6;
    elseif (lambda(1) == 0) && (lambda(1) == 0)
        indxConfig = 1;
    else
        indxConfig = 6 + 2*(lambda(1)-1);
    end
    if (alreadyDone(indxConfig) == 0)
        alreadyDone(indxConfig) = 1;
    else
        oldLambda2 = lambda;
        % find next available config
        if nbIter == nbIterMax
            break
        end
        jj = mod(indxConfig, nbIterMax) +1;
        while 1
            if (alreadyDone(jj) == 0)
                alreadyDone(jj) = 1;
                break
            else
                jj = mod(jj, nbIterMax) +1;
            end
        end
        
        indxZero = 0;
        switch(jj)
            case 1
                lambda = [0 0];
            case 2
                lambda = [-1 0];
                indxZero = 2;
            case 3
                lambda = [-1 -1];
            case 4
                lambda = [0 -1];
                indxZero = 1;
            case 5
                lambda = [1 -1];
            case 6 
                lambda = [1 0];
                indxZero = 2;
            case 7
                lambda = [-1 1];
            case 8
                lambda = [0 1];
                indxZero = 1;
            case 9
                lambda = [1 1];
        end
    end
end

if (nbIter == nbIterMax)
    keyboard;
end
end