function [atoms,bonds] = GenerateAtoms()
% function atoms = GenerateAtoms()
%
% Generates a list of bonds and atoms which correspond to a cyclohexane
% molecule near its equilibrium position. It is slightly out of equilibrium
% to see vibrations in the molecule which is relevant to the paper this
% project is based on.

zVal = 1.11;
atoms(18,1) = Atom;
bonds(18,1) = Bond;
for i = 0:5
    angle = 2*pi()/6;
    z = .4*((-1)^(mod(i, 2)));
    % Define atoms
    atoms(3*i+1).pos(1) = 1.5*cos(angle*i);
    atoms(3*i+1).pos(2) = 1.5*sin(angle*i);
    atoms(3*i+1).pos(3) = z;
    atoms(3*i+1).weight = 12;
    atoms(3*i+1).id = 3*i+1;
    atoms(3*i+1).force = [0,0,0];
    if z < 0 
        atoms(3*i+2).pos(1) = 2.25*cos(angle*i);
        atoms(3*i+2).pos(2) = 2.25*sin(angle*i);
        atoms(3*i+2).pos(3) = - z;
        atoms(3*i+2).weight = 1;
        atoms(3*i+2).id = 3*i+2;
        atoms(3*i+2).force = [0,0,0];
        
        atoms(3*i+3).pos(1) = 1.5*cos(angle*i);
        atoms(3*i+3).pos(2) = 1.5*sin(angle*i);
        atoms(3*i+3).pos(3) = z-zVal;
        atoms(3*i+3).weight = 1;
        atoms(3*i+3).id = 3*i+3;
        atoms(3*i+3).force = [0,0,0];
    else
        atoms(3*i+2).pos(1) = 1.5*cos(angle*i);
        atoms(3*i+2).pos(2) = 1.5*sin(angle*i);
        atoms(3*i+2).pos(3) = z+zVal;
        atoms(3*i+2).weight = 1;
        atoms(3*i+2).id = 3*i+2;
        atoms(3*i+2).force = [0,0,0];
        
        atoms(3*i+3).pos(1) = 2.25*cos(angle*i);
        atoms(3*i+3).pos(2) = 2.25*sin(angle*i);
        atoms(3*i+3).pos(3) = -z;
        atoms(3*i+3).weight = 1;
        atoms(3*i+3).id = 3*i+3;
        atoms(3*i+3).force = [0,0,0];
    end
    % Define bonds between C and H's
    bonds(3*i+1).atom1 = atoms(3*i+1);
    bonds(3*i+1).atom2 = atoms(3*i+2);
    bonds(3*i+2).atom1 = atoms(3*i+1);
    bonds(3*i+2).atom2 = atoms(3*i+3);
end
for i = 0:5
    % Define bonds between C's.
    bonds(3*i+3).atom1 = atoms(3*i+1);
    bonds(3*i+3).atom2 = atoms(mod(3*(i+1)+1,18));
end
    
    