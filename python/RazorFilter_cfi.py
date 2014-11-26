import FWCore.ParameterSet.Config as cms

RazorFilter = cms.EDFilter('RazorFilter' ,
                            jetInputTag = cms.InputTag("ak5PFJets"),
                            metInputTag = cms.InputTag("pfMet"),
                            btagInputTag = cms.InputTag("combinedSecondaryVertexMVABJetTags"),
                            csvFileName = cms.string("razor.csv"),
                            rootFileName = cms.string("razor.root"),
                            minJetPt = cms.double(30.0),
                            maxJetEta = cms.double(3.0)
                            )
