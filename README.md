# MuonHLTAnalyser
An analyser for the muon objects created by the HLT

This code should compile out of the box. Simply take a HLT menu and add the following at the end:
process.MuonHLTAna = cms.EDAnalyzer('MuonHLTAnalyser')
process.go = cms.Path(process.MuonHLTAna)
process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
