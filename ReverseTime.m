function atoms = ReverseTime(atoms, bonds, time)
% This function generates the atomic positions at the previous step in
% time.
% 
% atoms is a vector of atoms and time is the length of the time step in
% seconds

numAtoms = size(atoms,1);
for i = 1:numAtoms
    atoms(i).posBefore = atoms(i).pos - atoms(i).vel *...
        time + 0.5 * atoms(i).force/(atoms(i).weight * 10^-22)*...
        time^2;
end
end
