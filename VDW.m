function atoms = VDW(atoms, bonds)
% function atoms = VDW(atoms, bonds)
%
% This function calculates the van der walls forces for each atom and adds
% this to the force already experienced by each atom.

global PE;

numAtoms = size(atoms,1);
for i = 1:numAtoms - 1
    for j = i + 1:numAtoms
        numBonds = size(bonds,1);
        isBond = false;
        for k = 1:numBonds
            a = atoms(i).id;
            b = atoms(j).id;
            c = bonds(k).atom1.id;
            d = bonds(k).atom2.id;
            if a == c && b == d || ...
                    a == d && b == c
                isBond = true;
                break;
            end
        end
        if ~isBond    
            if atoms(i).weight == 12
                eps1 = .027;
                r1 = 2.04;
            else 
                eps1 = .020 ;
                r1 = 1.62;
            end
            if atoms(j).weight == 12
                eps2 = .027 ;
                r2 = 2.04;
            else
                eps2 = .020;
                r2 = 1.62;
            end
            eps = (eps1 + eps2)/2;
            r = (r1 + r2)/2;
            dist = atoms(j).pos - atoms(i).pos;
            dir = dist / norm(dist);
            d = norm(dist);
            force = eps * (13.5 * (r^6/d^7) - 22.08 * 10^5/r * ...
                exp(-12 * d/r));
            atoms(i).force = atoms(i).force + force*dir;
            atoms(j).force = atoms(j).force + force*dir;
            
            % Calculate potential energy
            PE = PE + eps * (-2.25 * (r/d)^6 + 1.84 * 10^5 *...
                exp(-12 * d/r)); % in kcal/mol
        end
    end
end
end
        
        
        