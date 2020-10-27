function atoms = Torsion(atoms, bonds)
% function atoms = AngleBending(atoms, bonds)
% 
% This function calculates the forces due to angle bending on all atoms and
% updates their forces. Note that theta0 is type 3 for all angles
% considered in a ring.

global PE;

numBonds = size(bonds,1);
numAngles = 0;
% Construct vector of angles which each contain 3 atoms.
for i = 1:numBonds - 1
    for j = i+1:numBonds
        atom1 = bonds(i).atom1;
        atom2 = bonds(i).atom2;
        atom3 = bonds(j).atom1;
        atom4 = bonds(j).atom2;
        if atom1.id == atom3.id || atom1.id == atom4.id || ...
                atom2.id == atom3.id || atom2.id == atom4.id
            numAngles = numAngles+1;
            angles(numAngles,:) = [bonds(i), bonds(j)];
        end
    end
end
numTorsions = 0;
% loop over angles to determine torsion groups.
for i = 1:numAngles-1
    for j = i+1
        if angles(i,1) == angles(j,1)
            numTorsions = numTorsions + 1;
            torsion(numTorsions,:) = [angles(i,2),angles(i,1),angles(j,2)];
        elseif angles(i,2) == angles(j,1)
            numTorsions = numTorsions + 1;
            torsion(numTorsions,:) = [angles(i,1),angles(i,2),angles(j,2)];
        elseif angles(i,1) == angles(j,2)
            numTorsions = numTorsions + 1;
            torsion(numTorsions,:) = [angles(i,2),angles(i,1),angles(j,1)];
        elseif angles(i,2) == angles(j,2)
            numTorsions = numTorsions + 1;
            torsion(numTorsions,:) = [angles(i,1),angles(i,2),angles(j,1)];
        end
    end
end
% loop over torsions to calculate forces
for i = 1:numTorsions
    if torsion(i,1).atom1 == torsion(i,2).atom1
        atom1 = torsion(i,1).atom2;
        atom2 = torsion(i,1).atom1;
        if torsion(i,2).atom2 == torsion(i,3).atom1
            atom3 = torsion(i,3).atom1;
            atom4 = torsion(i,3).atom2;
        else
            atom3 = torsion(i,3).atom2;
            atom4 = torsion(i,3).atom1;
        end
    elseif torsion(i,1).atom1 == torsion(i,2).atom2
        atom1 = torsion(i,1).atom2;
        atom2 = torsion(i,2).atom2;
        if torsion(i,2).atom1 == torsion(i,3).atom1
            atom3 = torsion(i,3).atom1;
            atom4 = torsion(i,3).atom2;
        else
            atom3 = torsion(i,3).atom2;
            atom4 = torsion(i,3).atom1;
        end
    elseif torsion(i,1).atom2 == torsion(i,2).atom1
        atom1 = torsion(i,1).atom1;
        atom2 = torsion(i,1).atom2;
        if torsion(i,2).atom2 == torsion(i,3).atom1
            atom3 = torsion(i,3).atom1;
            atom4 = torsion(i,3).atom2;
        else
            atom3 = torsion(i,3).atom2;
            atom4 = torsion(i,3).atom1;
        end
    else
        atom1 = torsion(i,1).atom1;
        atom2 = torsion(i,1).atom2;
        if torsion(i,2).atom1 == torsion(i,3).atom1
            atom3 = torsion(i,3).atom1;
            atom4 = torsion(i,3).atom2;
        else
            atom3 = torsion(i,3).atom2;
            atom4 = torsion(i,3).atom1;
        end
    end
    if size(unique([atom1, atom2, atom3, atom4]),2) == 4
        if atom1.weight + atom2.weight == 24
            % C-C-C-C 
            v1 = .185;
            v2 = .170;
            v3 = .520;
        elseif atom1.weight + atom2.weight == 13
            % C-C-C-H
            v1 = 0;
            v2 = 0;
            v3 = .280;
        else
            % H-C-C-H
            v1 = 0;
            v2 = 0;
            v3 = .238;
        end
        % Calculate torsion angle.
        a = atom2.pos - atom1.pos;
        b = atom3.pos - atom2.pos;
        c = atom4.pos - atom3.pos;
        n1 = cross(a,b);
        n2 = cross(b,c);
        angle = acosd(dot(n1,n2)/(norm(n1)*norm(n2)));
        angle = real(angle);
        
        force = (-v1/2*sind(angle)+v2*sind(angle)-3*v3/2*sind(3*angle));
        
        PE = PE + (v1/2) * (1 + cosd(angle)) +...
            (v2/2) * (1 - cosd(2 * angle)) +...
            (v3/2) * (1 + cosd(3 * angle)); % in kcal/mol
        if n1 ~= 0
            n1 = n1/norm(n1);
        end
        if n2 ~= 0
            n2 = n2/norm(n2);
        end
        atom1.force = atom1.force + n1*force;
        atom4.force = atom4.force + n2*force;
    end
end

