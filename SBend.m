function atoms = SBend(atoms, bonds)
% function atoms = SBend(atoms, bonds)
% 
% This function calculates the stretch-bend force on each atom
%
% atoms is a vector of atoms and bonds is a vector of bonds.

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
for i = 1:numAngles
    % determine type of angle
    if angles(i,1).weight + angles(i,2).weight + angles(i,3).weight == 36
        %C-C-C angle
        k = .13;
        theta0 = 111.0;
        l1 = 1.5247;
        l2 = 1.5247;
    elseif angles(i,1).weight + angles(i,3).weight == 2
        % H-C-H angle
        k = .08;
        theta0 = 109.5;
        l1 = 1.112;
        l2 = 1.112;
    else
        % C-C-H angle
        k = 0.;
        theta0 = 110.7;
        if angles(i,1).weight == 12
            l1 = 1.5247;
            l2 = 1.112;
        else
            l1 = 1.112;
            l2 = 1.5247;
        end
    end
    a = angles(i,1).pos - angles(i,2).pos;
    b = angles(i,3).pos - angles(i,2).pos;
    angle = acosd(dot(a,b)/(norm(a)*norm(b)));
    
    
    % Calculate stretch force
    force = 2.51118*k*(angle-theta0);
    dir = angles(i,3).pos - angles(i,2).pos;
    dir = dir/norm(dir);
    angles(i,3).force = angles(i,3).force - force * dir;
    angles(i,2).force = angles(i,2).force + force * dir;
    dir = angles(i,1).pos - angles(i,2).pos;
    dir = dir/norm(dir);
    angles(i,1).force = angles(i,1).force - force * dir;
    angles(i,2).force = angles(i,2).force + force * dir;
    
    % Calculate potential energy in kcal/mol
    PE = PE + 2.51118 * k * (((abs(norm(a)) - l1) +...
        (abs(norm(b)) - l2)) * (angle - theta0));
end
end