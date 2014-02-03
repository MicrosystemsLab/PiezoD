import sys, string, subprocess

# Parameter values
dopantValues = ['arsenic', 'boron', 'phosphorus']
doseValues = [2e14, 2e15, 2e16]
energyValues = [20, 50, 80]
timeValues = [15, 30, 45, 60, 75, 90, 105, 120]
tempValues = [900, 1000, 1100]

# Clear out all of the old simulation files
subprocess.call('rm -rf *.out', shell=True)
subprocess.call('rm -rf simulationInputs/*', shell=True)
subprocess.call('rm -rf simulationOutputs/*', shell=True)

# Generate the tsuprem input files and run tsuprem on them
ii = 1
for dopantLong in dopantValues:
  if dopantLong == 'boron':
    dopantShort = 'B'
    dopantBackground = 'phosphorus'
  elif dopantLong == 'phosphorus':
    dopantShort = 'P'
    dopantBackground = 'boron'
  elif dopantLong == 'arsenic':
    dopantShort = 'As'
    dopantBackground = 'boron'
  else:
    print 'Unknown dopant type!'
    sys.exit(0)

  for dose in doseValues:
    for energy in energyValues:
      for time in timeValues:
        for temp in tempValues:

          print 'Iteration %d of %d' % (ii, 2*len(dopantValues)*
            len(doseValues)*len(energyValues)*len(timeValues)*len(tempValues))
          ii += 1

          # First simulate without passivation oxide
          dryOxidationTime = 0
          wetOxidationTime = 0
          oxideStripCmd = ''

          # Build the simulation input file through substitution
          fileName = '%s_%s_%s_%s_%s_noOxide' % (dopantShort,
						dose, energy, temp, time)
          templateString = open('simulation.template', 'r').read()
          s = string.Template(templateString)
          fileContents = s.substitute(outputFileName=fileName, dose=dose,
            energy=energy, time=time, temp=temp, dopantShort=dopantShort,
            dopantLong=dopantLong, dopantBackground=dopantBackground,
            dryOxidationTime=dryOxidationTime, 
						wetOxidationTime=wetOxidationTime,
            oxideStripCmd = oxideStripCmd)

          # Save the file to disk
          simulationFileName = ''.join(['simulationInputs/', fileName, '.input'])
          f = open(simulationFileName, 'w')
          f.write(fileContents)
          f.close()

          # Run the simulation
          cmd = ''.join(['tsuprem4 ', simulationFileName])
          subprocess.call(cmd, shell=True)


          # Resimulate with an oxidation step immediately before the anneal
          dryOxidationTime = 5
          if temp == 900:
            wetOxidationTime = 66
          elif temp == 1000:
            wetOxidationTime = 15
          elif temp == 1100:
            wetOxidationTime = 5
          oxideStripCmd = 'etch oxide all'

          # Build the simulation input file through substitution
          fileName = '%s_%s_%s_%s_%s_oxide' % (dopantShort, dose,
						energy, temp, time)
          templateString = open('simulation.template', 'r').read()
          s = string.Template(templateString)
          fileContents = s.substitute(outputFileName=fileName, dose=dose,
            energy=energy, time=time, temp=temp, dopantShort=dopantShort,
            dopantLong=dopantLong, dopantBackground=dopantBackground,
            dryOxidationTime=dryOxidationTime, wetOxidationTime=wetOxidationTime,
            oxideStripCmd = oxideStripCmd)

          # Save the file to disk
          simulationFileName = ''.join(['simulationInputs/', fileName, '.input'])
          f = open(simulationFileName, 'w')
          f.write(fileContents)
          f.close()

          # Run the simulation
          cmd = ''.join(['tsuprem4 ', simulationFileName])
          subprocess.call(cmd, shell=True)

