# MuonHLTAnalyser
An analyser for the muon objects created by the HLT

This code should compile out of the box. Currently it is set to look at objects passing the Mu40 path. Changing the filter names would allow you to look at different paths.

Simply take a HLT menu (better to specifiy only the path you are interested in, rather than the whole menu):
hltGetConfiguration /dev/CMSSW_7_2_1/HLT/V56\
 --full --offline --mc --timing --unprescale --process BRSHLT\
 --globaltag auto:startup_GRun\
 --paths "HLTriggerFirstPath,HLT_Mu40_v*,HLTriggerFinalPath"\
 --output minimal --max-events -1 > hlt_Mu40_721.py

Then add the following at the end of the menu:
process.MuonHLTAna = cms.EDAnalyzer('MuonHLTAnalyser')
process.go = cms.Path(process.MuonHLTAna)
process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
