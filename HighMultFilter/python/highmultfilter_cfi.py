import FWCore.ParameterSet.Config as cms
highMultFilter = cms.EDFilter('HighMultFilter',
                                  multMax = cms.double(100000),
                                  multMin = cms.double(120),
                                  etaMax = cms.double(2.4),
                                  etaMin = cms.double(-2.4),
                                  doGenParticle = cms.bool(False),
                                  genSrc = cms.InputTag("genParticles"),
                                  vertexSrc = cms.InputTag("offlinePrimaryVertices"),
                                  trackSrc = cms.InputTag("generalTracks")
                                  )
