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
options.register('trackCollection',
                 'muonTracks',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Track collection: muonTracks, muonGlobalTracks, or gedElectronTracks")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddra_validation.root')
options.parseArguments()

# Map parallelRun.sh --track-collection names to HyddraLeptonTrackProducer config
_COLLECTION_MAP = {
    'muonTracks':       ('muon',     'slimmedMuons'),
    'muonGlobalTracks': ('muon',     'slimmedMuons'),
    'gedElectronTracks':('electron', 'slimmedElectrons'),
}
if options.trackCollection not in _COLLECTION_MAP:
    raise RuntimeError("Unknown trackCollection '%s'. Valid: %s"
                       % (options.trackCollection, list(_COLLECTION_MAP)))
_leptonType, _src = _COLLECTION_MAP[options.trackCollection]

process = cms.Process("HYDDRAVAL")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles or ['file:input.root']))

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

process.TFileService = cms.Service("TFileService",
    fileName=cms.string(options.outputFile))

# Track producer
process.load("RecoVertex.HyddraSVProducer.hyddraEXO_cfi")

process.hyddraLeptonTracks.leptonType = cms.string(_leptonType)
process.hyddraLeptonTracks.src        = cms.InputTag(_src)
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
    tracks              = cms.InputTag("hyddraLeptonTracks"),
    genParticles        = cms.InputTag("prunedGenParticles"),
    MET                 = cms.InputTag("slimmedMETs"),
    hasGenInfo          = cms.bool(options.hasGenInfo),
    motherPdgId         = cms.int32(options.motherPdgId),
    genDRCut            = cms.double(0.05),
    passSelDRCut        = cms.double(0.02),
)

process.p = cms.Path(process.hyddraLeptonTracks + process.hyddraSVsEXOProducer + process.hyddraVal)
process.schedule = cms.Schedule(process.p)
