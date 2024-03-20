from gecoris import dorisUtils
from gecoris import ioUtils
from pathlib import Path

outDir = '/home/caroline-pbazzocchi/algorithm/RC_ThesisTest_approx/RadarCoordinates/'

stationFiles = [f for f in Path('/home/caroline-pbazzocchi/algorithm/RC_ThesisTest_approx/').glob('????.json')]
stackFiles = [f for f in Path('/home/caroline-pbazzocchi/algorithm/RC_ThesisTest_approx/').glob('stack*.json')]


stations = []
for i in range(len(stationFiles)):
    stations.append(ioUtils.fromJSON(stationFiles[i]))
stacks = []
for i in range(len(stackFiles)):
    stacks.append(ioUtils.fromJSON(stackFiles[i]))

    
dorisUtils.RCexport(stations,stacks,outDir,plotFlag=2)
print('OK')
