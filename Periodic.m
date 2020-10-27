function periodicMatrix = Periodic(initialMatrix, b)
% periodicMatrix = Periodic(initialMatrix)
%
% initialMatrix: a matrix of integers (may be positive or negative)
%
% b: a positive integer indicating the upper bound requested
%
% periodicMatrix: a matrix of integers on the interval [1,b]
periodicMatrix = initialMatrix;
while min(min(periodicMatrix)) < 1 || max(max(periodicMatrix)) > b
    lowFilt = periodicMatrix < 1;
    periodicMatrix(lowFilt) = periodicMatrix(lowFilt) + b;
    highFilt = periodicMatrix > b;
    periodicMatrix(highFilt) = periodicMatrix(highFilt) - b;
end