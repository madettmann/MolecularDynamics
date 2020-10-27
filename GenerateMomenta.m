function [momenta] = GenerateMomenta(atoms, temperature)
% This function generates momenta for each of the atoms given in lattice.
%
% momenta is an array of momenta for each atom in lattice. Returned units
% are in gram-angstroms per second.
% 
% lattice: a list of atomic positions. Here, unit is arbitrary. 
% Temperature is the temperature of the system given in kelvin. 

% Uncomment for zero starting momenta
% numAtoms = size(lattice,1);
% momenta = zeros(numAtoms,3);


numAtoms = size(atoms, 1);
kB = 1.3806 * 10 ^(-23); % Boltzmann's constant m^2 . kg . s^(-2) . K^(-1)
% generate 3 random numbers sampled from a Gaussian distribution
randX = randn(numAtoms, 1);
randY = randn(numAtoms, 1);
randZ = randn(numAtoms, 1);

% find the average momentum of each atom and subtract from momenta of each
averageMomentumX = sum(randX, 1)/numAtoms;
averageMomentumY = sum(randY, 1)/numAtoms;
averageMomentumZ = sum(randZ, 1)/numAtoms;
randX = randX - averageMomentumX;
randY = randY - averageMomentumY;
randZ = randZ - averageMomentumZ;

sumMomenta = sum((randX.^2 + randY.^2 + randZ.^2), 1);
sumMass = 0;
for i = 1:numAtoms
    sumMass = sumMass + atoms(i).weight * 10^-25;
end

% calculate scaling factor, alpha.
alpha = sqrt((3/2*numAtoms*kB*temperature)./(1/2*sumMomenta/...
  (sumMass/numAtoms))); % kg-m/s
momenta = zeros(numAtoms, 3);
momenta(:, 1) = alpha .* randX;
momenta(:, 2) = alpha .* randY;
momenta(:, 3) = alpha .* randZ;

momenta = momenta * 10^13; % convert momenta to g-angstroms per second.

end
