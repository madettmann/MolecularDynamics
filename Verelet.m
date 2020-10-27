function positionsAfter = Verelet(atoms, bonds, timeStep)
% This functions implements the Verelet algorithm to calculate a list of
% positions at time t + dt.

% atoms is a vector of atoms, bonds is a vector of bonds, timeStep is the
% length of the time step in seconds. it returns positionsAfter which is a
% matrix containing the new positions in the form (x y z)

% calculate the positions at t+dt
numAtoms = size(atoms,1);
positionsAfter = zeros(numAtoms,3);
for i = 1:numAtoms
    positionsAfter(i,:) = 2 .* atoms(i).pos - atoms(i).posBefore +...
        atoms(i).force/(atoms(i).weight*10^-22) .* timeStep ^2;
end
end