#ifndef MuonHLTAnalyser_H
#define MuonHLTAnalyser_H

// Default
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Additional
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "Math/GenVector/VectorUtil.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"

// Physics Objects
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

//Trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"


// ROOT
#include "TH1.h"
#include "TMath.h"

#include <string>
#include <fstream>

class MuonHLTAnalyser : public edm::EDAnalyzer {
public:
	explicit MuonHLTAnalyser(const edm::ParameterSet&);
	~MuonHLTAnalyser();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

//	Num of Seeds:
	TH1D *histNumL3IOHitSeeds;

	TH1D *histhltL3TrajSeedOIState;
	TH1D *histhltL3TrajSeedOIHit;
	TH1D *histhltL3TrajSeedIOHit;

//	Num of Tracker Tracks:
	TH1D *histNumL3TkTrkFromL2IOHit;
	TH1D *histNumTkTrkIOIter;
	TH1D *histhltL3TkTracksFromL2OIState;
	TH1D *histhltL3TkTracksFromL2OIHit;
	TH1D *histhltL3TkTracksFromL2IOHit;


//	Eff Distributions:
	TH1D *hGenEta;
	TH1D *hGenPhi;
	TH1D *hGenPt;

	TH1D *hAllL1Eta;
	TH1D *hAllL1Phi;
	TH1D *hAllL1Pt;

	TH1D *hAllL2Eta;
	TH1D *hAllL2Phi;
	TH1D *hAllL2Pt;
	TH1D *hPassL2Eta;
	TH1D *hPassL2Phi;
	TH1D *hPassL2Pt;

	TH1D *hNumL1MatchedToGen;
	TH1D *hNumL2MatchedToGen;
	TH1D *hNumL3MatchedToGen;
	TH1D *hNumAllL3MatchedToGen;
	TH1D *hNumL1MatchedToGenInEvent;
	TH1D *hNumL2MatchedToGenInEvent;
	TH1D *hNumL3MatchedToGenInEvent;
	TH1D *hNumAllL3MatchedToGenInEvent;

	TH1D *hL1Eta;
	TH1D *hL1Pt;
	TH1D *hL2Eta;
	TH1D *hL2Pt;

	TH1D *hNumL2Found;
	TH1D *hNumL3Found;
	TH1D *hNumNoL3Found;

	TH1D *hL1DistToGen;
	TH1D *hL2DistToGen;
	TH1D *hL3DistToGen;

	TH1D *hNumGenPerEvent;

	TH1D *hAllDoubleL2Eta;
	TH1D *hPassDoubleL2Eta;
	TH1D *hAllDoubleL1Eta;
	TH1D *hPassDoubleL1Eta;

	TH1D *hAllDoubleL2Pt;
	TH1D *hPassDoubleL2Pt;
	TH1D *hAllDoubleL1Pt;
	TH1D *hPassDoubleL1Pt;

	TH1D *hL1DrGen;
	TH1D *hL2DrGen;
	TH1D *hL3DrGen;

	unsigned int numEvents;
	edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l2Token_;
};
#endif
