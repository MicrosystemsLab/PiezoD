function postProcessTables()
clear all
close all
clc

% size = (3, 3, 3, 3, 8, 2) -> (dopant, dose, energy, temp, time, oxidation)
dopants = {'B', 'P', 'As'};
dopantsLong = {'boron', 'phosphorus', 'arsenic'};
doses = [2e14, 2e15, 2e16];
energies = [20 50 80];
temps = [900 1000 1100];
times = [15 30 45 60 75 90 105 120];

% Setup the z interpolation grid
% All data will be calculated over this grid
z = 0:0.01:5;

% Allocate space
Beta1 = zeros(3, 3, 3, 3, 8, 2);
Beta2 = zeros(3, 3, 3, 3, 8, 2);
Rs = zeros(3, 3, 3, 3, 8, 2);
Nz_total = zeros(3, 3, 3, 3, 8, 2);
Nz = zeros(3, 3, 3, 3, 8, 2);
Xj = zeros(3, 3, 3, 3, 8, 2);
n = zeros(length(z), 3, 3, 3, 3, 8, 2);

ImplantDopants = [1 2 3];
ImplantEnergies = energies;
ImplantDoses = doses;
AnnealTemps = temps;
AnnealTimes = times;
AnnealOxidation = [1 2];

currentCount = 1;
for ii = 1:length(dopants)
  for jj = 1:length(doses)
    for kk = 1:length(energies)
      for ll = 1:length(temps)
        for mm = 1:length(times)
          dopant = dopants{ii};
          dopantLong = dopantsLong{ii};
          dose = doses(jj);
          energy = energies(kk);
          temp = temps(ll);
          time = times(mm);
          
          % Output a status update
          fprintf('Processing %s, %0.1g/sq cm, %dkeV, %dC, %dmin\n', ...
						dopant, dose, energy, temp, time);
          fprintf('%d of %d\n\n', currentCount, 2*length(dopants)* ...
						length(doses)*length(energies)*length(temps)*length(times));
          
          % Generate the filenames
          doseExponent = floor(log10(dose));
          dosePrefactor = round(10^(rem(log10(dose), doseExponent)));
          name = sprintf('%s_%de+%d_%d_%d_%d', dopant, ...
						dosePrefactor, doseExponent, energy, temp, time);
          filename_noOxide = ['simulationOutputs/' name '_noOxide.out'];
          filename_oxide = ['simulationOutputs/' name '_oxide.out'];

					% No oxide case
          fprintf('%s\n', filename_noOxide);
          nInterp = loadAndCleanupData(filename_noOxide);
          [betatmp, Beta1tmp, Beta2tmp, Nztmp, Nz_totaltmp, Rstmp, Xjtmp] = ...
						calculateProfileProperties(dopantLong, z*1e-6, nInterp);
          Rs(ii, jj, kk, ll, mm, 1) = Rstmp;
          Xj(ii, jj, kk, ll, mm, 1) = Xjtmp;
          Beta1(ii, jj, kk, ll, mm, 1) = Beta1tmp;
          Beta2(ii, jj, kk, ll, mm, 1) = Beta2tmp;
          Nz(ii, jj, kk, ll, mm, 1) = Nztmp;
          Nz_total(ii, jj, kk, ll, mm, 1) = Nz_totaltmp;
          n(:, ii, jj, kk, ll, mm, 1) = nInterp;

					% Oxide case
          fprintf('%s\n', filename_oxide);
          nInterp = loadAndCleanupData(filename_oxide);
          [betatmp, Beta1tmp, Beta2tmp, Nztmp, Nz_totaltmp, Rstmp, Xjtmp] = ...
						calculateProfileProperties(dopantLong, z*1e-6, nInterp);
          Rs(ii, jj, kk, ll, mm, 2) = Rstmp;
          Xj(ii, jj, kk, ll, mm, 2) = Xjtmp;
          Beta1(ii, jj, kk, ll, mm, 2) = Beta1tmp;
          Beta2(ii, jj, kk, ll, mm, 2) = Beta2tmp;
          Nz(ii, jj, kk, ll, mm, 2) = Nztmp;
          Nz_total(ii, jj, kk, ll, mm, 2) = Nz_totaltmp;
          n(:, ii, jj, kk, ll, mm, 2) = nInterp;
          
          % Move to the next iteration
          currentCount = currentCount + 2;
        end
      end
    end
  end
end

ratio = Nz./Nz_total;
ratio = reshape(ratio, 1, 3*3*3*3*8*2);
fprintf('Nz/Nz_total: %.2f to %.2f (mean/median/std = %.2f/%.2f/%.2f)\n', ...
	min(ratio), max(ratio), mean(ratio), median(ratio), std(ratio));

% Output data
save lookupTable z Rs Xj Beta1 Beta2 Nz Nz_total n ImplantDopants ...
  ImplantEnergies ImplantDoses AnnealTemps AnnealTimes AnnealOxidation;

  function nInterp = loadAndCleanupData(filename)
    [zRaw, nRaw, material] = textread(filename, '%f %f %s', 'headerlines', 1);
    siliconIndices = find(strcmp(material, 'silicon'));
    zRaw = zRaw(siliconIndices);
    nRaw = nRaw(siliconIndices);
    zRaw = zRaw - zRaw(1);
    nInterp = interp1(zRaw, nRaw, z, 'pchip'); % Arrange data on a standard grid
  end
end