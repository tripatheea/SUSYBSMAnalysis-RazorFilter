import FWCore.ParameterSet.Config as cms

DijetFilter = cms.EDFilter("SimpleJetFilter",
                         jetCollection = cms.InputTag("ak5PFJets"),
                         ptCut = cms.double(60.),
                         maxRapidityCut = cms.double(3.0),
                         nJetMin = cms.uint32(2)
                         )

