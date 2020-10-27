function atoms = BBend(atoms, bonds)
% function atoms = BBend(atoms, bonds)
%
% This function calculates the bend-bend forces on each atom.

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
        if atom1.id == atom3.id
            numAngles = numAngles+1;
            angles(numAngles,:) = [atom2, atom1, atom4];
        elseif atom1.id == atom4.id
            numAngles = numAngles+1;
            angles(numAngles,:) = [atom2, atom1, atom3];
        elseif atom2.id == atom3.id
            numAngles = numAngles+1;
            angles(numAngles,:) = [atom1, atom3, atom4];
        elseif atom2.id == atom4.id
            numAngles = numAngles+1;
            angles(numAngles,:) = [atom1, atom2, atom3];
        end
    end
end
numBBend = 0;
for i = 1:numAngles-1
    for j = i+1:numAngles
        if angles(i,2).id == angles(j,2).id
            % Angles share center atom.
            if angles(i,1).weight + angles(i,3).weight == 24
                % C-C-C 
                k1 = .24;
                angle01 = 110.;
            elseif angles(i,1).weight + angles(i,3).weight == 13
                % C-C-H
                k1 = .30;
                angle01 = 110.7;
            else 
                k1 = 0.;
                angle01 = 109.5;
            end
            if angles(j,1).weight + angles(j,3).weight == 24
                % C-C-C 
                k2 = .24;
                angle02 = 110.;
            elseif angles(j,1).weight + angles(j,3).weight == 13
                % C-C-H
                k2 = .30;
                angle02 = 110.7;
            else 
                k2 = 0.;
                angle02 = 109.5;
            end
            k = k1*k2;
            a = angles(i,1).pos - angles(i,2).pos;
            b = angles(i,3).pos - angles(i,2).pos;
            angle1 = acosd(dot(a,b)/(norm(a)*norm(b)));
            n = cross(a,b); 
            dir1 = cross(n,a);
            if dir1 ~= 0
                dir1 = -dir1/norm(dir1);
            end
            a = angles(j,1).pos - angles(j,2).pos;
            b = angles(j,3).pos - angles(j,2).pos;
            angle2 = acosd(dot(a,b)/(norm(a)*norm(b)));
            dir2 = cross(n,a);
            if dir2 ~= 0
                dir2 = -dir2/norm(dir2);
            end
            force1 = -0.021914*k*(angle2 - angle02);
            force2 = -0.021914*k*(angle1 - angle01);
            
            % Calculate potential energy in kcal/mol
            PE = PE - .021914 * k *(angle2 - angle02) * (angle1 - angle01);
            
            angles(i,1).force = angles(i,1).force - force1*dir1;
            angles(i,3).force = angles(i,3).force + force1*dir1;
            
            angles(j,1).force = angles(j,1).force - force2*dir2;
            angles(j,3).force = angles(j,3).force + force2*dir2;
        end
    end
end
end