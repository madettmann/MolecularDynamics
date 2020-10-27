function atoms = AngleBending(atoms, bonds)
% function atoms = AngleBending(atoms, bonds)
% 
% This function calculates the forces due to angle bending on all atoms and
% updates their forces. Note that theta0 is type 3 for all angles
% considered in a ring.

global PE;

numBonds = size(bonds,1);
numAngles = 0;
% Construct vector of angles which each contain 3 atoms.
% Don't initialize angles vector outside loop because I wanted it to remain
% generic.
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
% loop over angles
for i = 1:numAngles
    % determine type of angle
    if angles(i,1).weight + angles(i,2).weight + angles(i,3).weight == 36
        %C-C-C angle
        k = .67;
        theta0 = 111.0;
    elseif angles(i,1).weight + angles(i,3).weight == 2
        % H-C-H angle
        k = .55;
        theta0 = 109.5;
    else
        % C-C-H angle
        k = .59;
        theta0 = 110.7;
    end
    a = angles(i,1).pos - angles(i,2).pos;
    b = angles(i,3).pos - angles(i,2).pos;
    angle = acosd(dot(a,b)/(norm(a)*norm(b)));
    force = -(.043828*k*(angle - theta0) * (1 - .014 * (angle-theta0) + ...
        5.6 * 10^-5 * (angle - theta0)^2 - 7 * 10^-7 * (angle-theta0)^3 ...
        + 9. * 10^-10 * (angle-theta0)^4) + .021914 * k * (angle-theta0)...
        * (.014+11.2 * 10^-5 * (angle-theta0) - 21 * 10^-7 * ...
        (angle-theta0)^2 + 36 * 10^-10 * (angle-theta0)^3));
    % Calculate potential energy in kcal/mol
    PE = PE + .021914 * k * (angle - theta0)^2 * (1 - .014 *...
        (angle - theta0) + 5.6 * 10^-5 * (angle - theta0)^2 -...
        7.0 * 10^-7 * (angle - theta0)^3 + 9.0 * 10^-10 * (angle - theta0));
    % define vector normal to angle.
    n = cross(a,b);
    dir = cross(n,a);
    if dir ~= 0
        dir = -dir/norm(dir);
    end
    angles(i,1).force = angles(i,1).force + force*dir;
    dir = cross(n,b) ;
    if dir ~= 0
        dir = -dir/norm(dir);
    end
    angles(i,3).force = angles(i,3).force - force*dir; % mdyne
end
end
