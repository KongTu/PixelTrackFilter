// -*- C++ -*-
//
// Package:    PixelTrackFilter/PixelTrackFilter
// Class:      PixelTrackFilter
// 
/**\class PixelTrackFilter PixelTrackFilter.cc PixelTrackFilter/PixelTrackFilter/plugins/PixelTrackFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Sun, 15 May 2016 13:11:33 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PixelTrackFilter : public edm::EDFilter {
public:
    explicit PixelTrackFilter(const edm::ParameterSet&);
    ~PixelTrackFilter();
    virtual void endJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
private:

    edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandSrc_;
    edm::EDGetTokenT<CaloTowerCollection> towerSrc_;
    
    double multMax_;
    double multMin_;
    double etaMax_;
    double etaMin_;

    bool doGenParticle_;
    bool doDS_;
    bool doDS_caloTower_;

    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PixelTrackFilter::PixelTrackFilter(const edm::ParameterSet& iConfig) :
genSrc_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
vertexSrc_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
trackSrc_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("trackSrc"))),
pfCandSrc_(consumes<reco::PFCandidateCollection >(iConfig.getParameter<edm::InputTag>("pfCandSrc"))),
towerSrc_(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("towerSrc"))),
multMax_(iConfig.getParameter<double>("multMax")),
multMin_(iConfig.getParameter<double>("multMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
etaMin_(iConfig.getParameter<double>("etaMin")),
doGenParticle_(iConfig.getParameter<bool>("doGenParticle")),
doDS_(iConfig.getParameter<bool>("doDS")),
doDS_caloTower_(iConfig.getParameter<bool>("doDS_caloTower"))
{
    
}


PixelTrackFilter::~PixelTrackFilter()
{

}

bool
PixelTrackFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    bool accepted = false;

    double nMult_ass_good = 0;

    if( doGenParticle_ ){

    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(genSrc_,genpars);
    
    for(unsigned it=0; it<genpars->size(); ++it){
        
        const reco::GenParticle & trk = (*genpars)[it];
        
        double eta = trk.eta();
        if(eta>etaMax_ || eta<etaMin_) continue;
        double pt  = trk.pt();
        if(pt<=0.4) continue;
        if(trk.status()!=1) continue;
        if(fabs(trk.charge())!=1) continue;
        
        nMult_ass_good++;
    }

    }
    else{

    Handle<edm::View<reco::Track>> tracks;
    iEvent.getByToken(trackSrc_, tracks);

      for(unsigned it = 0; it < tracks->size(); it++){

       const reco::Track & trk = (*tracks)[it];

          if( trk.algo() != 4 ) continue;
          if( trk.pt() > 0.4 ){nMult_ass_good++;}// NtrkOffline        

      } 
    }

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexSrc_,vertices);
    double bestvz=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); 
    bool validVertex = false;
    if( !vtx.isFake() && vtx.tracksSize() >= 2 && fabs(bestvz) < 15 ) validVertex = true; //valid vertex selection


    Handle<reco::PFCandidateCollection> pfCandidates;
    iEvent.getByToken(pfCandSrc_, pfCandidates);

    double towerPlus = 0.0;
    double towerMinus = 0.0;
    for( unsigned ic = 0; ic < pfCandidates->size(); ic++ ) {

        const reco::PFCandidate& cand = (*pfCandidates)[ic];
        double ecalEnergy = cand.ecalEnergy();
        double hcalEnergy = cand.hcalEnergy();

        if( ( ecalEnergy+hcalEnergy ) > 3.0 && cand.eta() > 3.0 && cand.eta() < 5.0 ) towerPlus++;
        if( ( ecalEnergy+hcalEnergy ) > 3.0 && cand.eta() > -5.0 && cand.eta() < -3.0 ) towerMinus++;

    }

    Handle<CaloTowerCollection> towers;
    iEvent.getByToken(towerSrc_, towers);

    
    double caloTowerPlus = 0.0;
    double caloTowerMinus = 0.0;
    for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];
        double ecalEnergy = cand.emEnergy();
        double hcalEnergy = cand.hadEnergy();
        
        if( ( ecalEnergy+hcalEnergy ) > 3.0 && hit.eta() > 3.0 && hit.eta() < 5.0 ) caloTowerPlus++;
        if( ( ecalEnergy+hcalEnergy ) > 3.0 && hit.eta() > -5.0 && hit.eta() < -3.0 ) caloTowerMinus++;

    }
    if( doDS_caloTower_ ){
        if( caloTowerPlus > 0.0 && caloTowerMinus > 0.0 ) accepted = true;
    }
    if( doDS_ ){
        if( towerPlus > 0.0 && towerMinus > 0.0 ) accepted = true;
    }
    else{
        if(nMult_ass_good>=multMin_ && nMult_ass_good<multMax_) accepted = true;
    }
    
    return accepted;
}
void
PixelTrackFilter::endJob()
{
}
DEFINE_FWK_MODULE(PixelTrackFilter);