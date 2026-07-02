#run folding code
from gwfolding import foldSID
import Fr
paramsFile = '/home/jishnu.suresh/O4_testbed/O4_refactor/foldSID.par'
foldedInvCov, foldedStatistic, params, foldParams = foldSID(paramsFile, writeFrames=1)
