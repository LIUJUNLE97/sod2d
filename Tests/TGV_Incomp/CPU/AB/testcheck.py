# Small script to check the TGV  test results

import numpy as np

# Define a tolerance for numerical comparison
tol = 1.0e-4

# Load the current analysis results in the .dat file
curRes = open('analysis_cube-4.dat', 'r')

# Load the reference results in the .dat file
refRes = open('Reference/analysis_cube-4.dat', 'r')

# Check that the files are matched in number of time-steps (number of line entries)
curLines = curRes.readlines()
refLines = refRes.readlines()

if len(curLines)!=len(refLines):
    print('Number of time-steps in current and reference results do not match')
    exit()

# For all 6 columns, check that the entries are similarly matched
for i in range(len(curLines)):
    curEntries = curLines[i].split()
    refEntries = refLines[i].split()
    for j in range(len(curEntries)):
        aux = abs(float(curEntries[j])-float(refEntries[j]))/(abs(float(refEntries[j]))+tol)
        if aux>tol:
            print(i, j, aux)
            print('Mismatch in current and reference results')
            exit()
            
# Close the reference and current results
curRes.close()
refRes.close()
