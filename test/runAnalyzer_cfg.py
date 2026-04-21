import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('hasGenInfo',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Process gen-level information (set False for data)")
options.register('motherPdgId',
                 9000006,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "PDG ID of the signal mother particle for gen matching")
options.register('maxNormChi2',
                 5.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Max vertex chi2/ndof")
options.register('outputCollection',
                 'seeds',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Base collection: seeds, inclusive, or isolated")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddra_validation.root')
options.parseArguments()

process = cms.Process("HYDDRAVAL")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles or ['file:input.root']))

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

process.TFileService = cms.Service("TFileService",
    fileName=cms.string(options.outputFile))

# Producer
process.load("RecoVertex.HyddraSVProducer.hyddraEXO_cfi")
process.hyddraSVsEXOProducer.leptonic.maxNormChi2 = cms.double(options.maxNormChi2)

# Analyzer
process.hyddraVal = cms.EDAnalyzer("HyddraSVsEXOAnalyzer",
    outputCollection    = cms.string(options.outputCollection),
    seedVertices        = cms.InputTag("hyddraSVsEXOProducer", "seedVertices"),
    inclusiveVertices   = cms.InputTag("hyddraSVsEXOProducer", "inclusiveVertices"),
    isolatedVertices    = cms.InputTag("hyddraSVsEXOProducer", "isolatedVertices"),
    disambiguationFlags = cms.InputTag("hyddraSVsEXOProducer", "disambiguationFlags"),
    seedIsolationFlags  = cms.InputTag("hyddraSVsEXOProducer", "seedIsolationFlags"),
    isolationFlags      = cms.InputTag("hyddraSVsEXOProducer", "isolationFlags"),
    pvCollection        = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks              = cms.InputTag("hyddraSVsEXOProducer", "leptonTracks"),
    genParticles        = cms.InputTag("prunedGenParticles"),
    MET                 = cms.InputTag("slimmedMETs"),
    hasGenInfo          = cms.bool(options.hasGenInfo),
    motherPdgId         = cms.int32(options.motherPdgId),
    genDRCut            = cms.double(0.05),
    passSelDRCut        = cms.double(0.02),
)

process.p = cms.Path(process.hyddraSVsEXOProducer + process.hyddraVal)
process.schedule = cms.Schedule(process.p)
