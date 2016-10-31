import FWCore.ParameterSet.Config as cms
pixelTrackFilter = cms.EDFilter('PixelTrackFilter',
                                  multMax = cms.double(100000),
                                  multMin = cms.double(1),
                                  etaMax = cms.double(2.4),
                                  etaMin = cms.double(-2.4),
                                  doGenParticle = cms.bool(False),
                                  doDS = cms.bool(False),
                                  doDS_caloTower = cms.bool(False),
                                  genSrc = cms.InputTag("genParticles"),
                                  vertexSrc = cms.InputTag("offlinePrimaryVertices"),
                                  trackSrc = cms.InputTag("generalTracks"),
                                  pfCandSrc = cms.InputTag("particleFlow"),
                                  towerSrc = cms.InputTag("towerMaker")
                                  )

