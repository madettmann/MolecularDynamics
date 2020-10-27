function atoms = Stretching(atoms, bonds)
% function atoms = CalculateForces(atoms, bonds)
%
% This function calculates stretching forces on each atom and adds them to
% the forces felt by each atom.

global PE;
numBonds = size(bonds,1);

for i = 1:numBonds
    weight1 = bonds(i).atom1.weight;
    weight2 = bonds(i).atom2.weight;
    pos1 = bonds(i).atom1.pos;
    pos2 = bonds(i).atom2.pos;
    l = norm(pos2-pos1);
    if weight1 + weight2 == 13 % C-H bond
        l0 = 1.112; 
        ks = 4.74; % Newtons/Angstroms
    else % C-C bond
        l0 = 1.5247;
        ks = 4.49; % Newtons/Angstroms
    end
    force = (143.88 * ks * (l - l0) * (1 - 2.55 * (l - l0) + 1.4875 * ...
        (l - l0) ^2) + 71.94 * ks * (l - l0) ^2 * (-2.55 * l + 2.98 *...
        (l - l0))); % in mdyne.
    % Calculate Potential energy
    PE = PE + 71.94 * ks * (l - l0)^2 * (1 - 2.55 * (l - l0) + 1.49 *...
        (l - l0)^2); % kcal/mol
    dir = (pos2 - pos1) / l;
    force1 = force * dir;
    bonds(i).atom1.force = bonds(i).atom1.force + force1;
    bonds(i).atom2.force = bonds(i).atom2.force - force1;
end