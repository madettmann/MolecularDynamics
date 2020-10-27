function ShowMolecule(atoms, bonds)
% function ShowMolecule(atoms, bonds)
%
% atoms is a vector of atoms and bonds is a vector of bonds.


atomLen = size(atoms,1);
bondLen = size(bonds,1);
x = zeros(atomLen,1);
y = zeros(atomLen,1);
z = zeros(atomLen,1);
c = zeros(atomLen,3);
% Draw atoms
for i = 1:atomLen
    x(i) = atoms(i).pos(1);
    y(i) = atoms(i).pos(2);
    z(i) = atoms(i).pos(3);
    % Color C atoms Red and H atoms blue
    if atoms(i).weight > 1
        c(i,:) = [1 0 0];
    else
        c(i,:) = [0 0 1];
    end
end
scatter3(x,y,z,1000,c,'filled');
% Add bonds between atoms.
for i = 1:bondLen
    pos1 = bonds(i).atom1.pos;
    pos2 = bonds(i).atom2.pos;
    line([pos1(1) pos2(1)],[pos1(2) pos2(2)], [pos1(3) pos2(3)],...
        'Color','black','LineWidth', 5);
end
xlim([-3 3]);
ylim([-3 3]);
zlim([-3 3]);
%axis equal