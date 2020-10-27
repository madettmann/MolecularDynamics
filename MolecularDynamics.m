% performs a molecular dynamics simulation and plots the  total energy as
% a function of time. This assumes that atoms interact with
% the MM3 potential.

steps = 10000; % number of steps in the simulation
ccLens = zeros(steps,1); % where the average C-C bond lengths are stored
chLens = zeros(steps,1); % where the average H-H bond lengths are stored
PEs = zeros(steps,1); % where the potential energies are stored
KEs = zeros(steps,1); % where the kinetic energies are stored


temperature = 0; % in Kelvin
dt = 10 ^-14; % in seconds
global PE;
KE = 0;
[atoms,bonds] = GenerateAtoms();

numAtoms = size(atoms, 1);
momenta = GenerateMomenta(atoms, temperature); % in g Angstroms/second
% Initialize velocities
for i = 1:numAtoms
    atoms(i).vel = momenta(i) / (atoms(i).weight*10^-22); % Angstroms per s
end

% Perform initialization of forces and calculate positions at -dt.
atoms = CalculateForces(atoms, bonds);
atoms = ReverseTime(atoms, bonds, dt);
positionsAfter = Verelet(atoms, bonds, dt);
% Update Velocities
for i = 1:numAtoms
    displacement = positionsAfter(i,:) - atoms(i).posBefore;
    atoms(i).vel = displacement / (dt * 2);
    KE = KE + .5 * atoms(i).weight * 10^-22 * (norm(atoms(i).vel))^2; % kcals
end

figure(1)
ShowMolecule(atoms,bonds);
f = waitbar(0,'Performing the Simulation');
numBonds = size(bonds,1);
for j = 1:numBonds
    % Calculate average CH bond lengths and CC bond lengths
  if bonds(j).atom1.weight + bonds(j).atom2.weight == 13
      chLens(1) = chLens(1) + norm(bonds(j).atom1.pos - ...
          bonds(j).atom2.pos);
  else
      ccLens(1) = ccLens(1) + norm(bonds(j).atom1.pos - ...
          bonds(j).atom2.pos);
  end
end
chLens(1) = chLens(1)/12;
ccLens(1) = ccLens(1)/6;
KEs(1) = KE;
PEs(1) = PE;
for i = 2:steps
    KE = 0;
    for j = 1:numAtoms
        atoms(j).posBefore = atoms(j).pos;
        atoms(j).pos = positionsAfter(j,:);
    end
    atoms = CalculateForces(atoms, bonds);
    positionsAfter = Verelet(atoms, bonds, dt);
    for j = 1:numAtoms
      displacement = positionsAfter(j,:) - atoms(j).posBefore;
      atoms(j).vel = displacement / (dt * 2);
      KE = KE + .5 * atoms(j).weight * 10^-22 * (norm(atoms(j).vel))^2;
    end
    if mod(i,100) == 0
        figure(1)
        ShowMolecule(atoms,bonds)
        pause(.0001);
        %saveas(gcf,strcat(num2str(i),'.png')); % Saves to a file for video
    end
    numBonds = size(bonds,1);
    for j = 1:numBonds
        if bonds(j).atom1.weight + bonds(j).atom2.weight == 13
            chLens(i) = chLens(i) + norm(bonds(j).atom1.pos - ...
                bonds(j).atom2.pos);
        else
            ccLens(i) = ccLens(i) + norm(bonds(j).atom1.pos - ...
                bonds(j).atom2.pos);
        end
    end
    
    chLens(i) = chLens(i)/12;
    ccLens(i) = ccLens(i)/6;
    KEs(i) = KE;
    PEs(i) = PE;
    waitbar(i/steps,f);
end
close(f);
figure(1)
ShowMolecule(atoms,bonds)
figure(2)
plot((1:steps)*dt,ccLens,(1:steps)*dt,chLens)
xlabel('Time (s)')
ylabel('Bond length (Angstroms)');
figure(3)
plot((1:steps)*dt,KEs + PEs);
xlabel('Time (s)')
ylabel('Total Energy (kcal/mol)');