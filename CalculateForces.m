function atoms = CalculateForces(atoms, bonds)
% function atoms = CalculateForces(atoms, bonds)
%
% This function, along with several helper functions, calculates the net
% force for each atom. There are several forces considered. Each is handled
% in turn.

global PE;

numAtoms = size(atoms,1);
for i = 1:numAtoms
    atoms(i).force = [0,0,0];
end

PE = 0;

% Calculate stretching forces
atoms = Stretching(atoms, bonds);
% Calculate angle bending forces
atoms = AngleBending(atoms, bonds);
% Calculate torsion forces
atoms = Torsion(atoms, bonds);
% Calculate bend-bend forces
atoms = BBend(atoms, bonds);
% Calculate van Der Waal forces
atoms = VDW(atoms, bonds);

% Calculate stretch-bend forces
atoms = SBend(atoms, bonds);



end