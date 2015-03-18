// Class:      MuonHLTAnalyser
// Description: Analyser for the Muon HLT performance
// Original Author:  Benjamin Radburn-Smith
// Created:  17 Mar 2015

#include "MuonHLTAnalyser.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

//HIT TEST
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"



MuonHLTAnalyser::MuonHLTAnalyser(const edm::ParameterSet& iConfig){
	l2Token_ = consumes<reco::RecoChargedCandidateCollection>(edm::InputTag("hltL2MuonCandidates","","BRSHLT"));
//	^ This is the proper way to get the collections ... but i still have the old ways later on

	edm::Service<TFileService> fs;

//	Num of Seeds:
	histNumL3IOHitSeeds = fs->make<TH1D>("histNumL3IOHitSeeds","histNumL3IOHitSeeds",50,0.0,50.0);

	histhltL3TrajSeedOIState = fs->make<TH1D>("histhltL3TrajSeedOIState","histhltL3TrajSeedOIState",50,0.0,50.0);
	histhltL3TrajSeedOIHit = fs->make<TH1D>("histhltL3TrajSeedOIHit","histhltL3TrajSeedOIHit",50,0.0,50.0);
	histhltL3TrajSeedIOHit = fs->make<TH1D>("histhltL3TrajSeedIOHit","histhltL3TrajSeedIOHit",50,0.0,50.0);


//	Num of Tracker Tracks:
	histNumL3TkTrkFromL2IOHit = fs->make<TH1D>("histNumL3TkTrkFromL2IOHit","histNumL3TkTrkFromL2IOHit",50,0.0,50.0);
	histNumTkTrkIOIter = fs->make<TH1D>("histNumTkTrkIOIter","histNumTkTrkIOIter",50,0.0,50.0);

	histhltL3TkTracksFromL2OIState = fs->make<TH1D>("histhltL3TkTracksFromL2OIState","histhltL3TkTracksFromL2OIState",50,0.0,50.0);
	histhltL3TkTracksFromL2OIHit = fs->make<TH1D>("histhltL3TkTracksFromL2OIHit","histhltL3TkTracksFromL2OIHit",50,0.0,50.0);
	histhltL3TkTracksFromL2IOHit = fs->make<TH1D>("histhltL3TkTracksFromL2IOHit","histhltL3TkTracksFromL2IOHit",50,0.0,50.0);

//	Eff Distributions:
	hGenEta = fs->make<TH1D>("hGenEta","hGenEta",100,-5.0,5.0);
	hGenPhi = fs->make<TH1D>("hGenPhi","hGenPhi",70,-3.5,3.5);
	hGenPt = fs->make<TH1D>("hGenPt","hGenPt",1000,0,1000.0);

	hAllL1Eta = fs->make<TH1D>("hAllL1Eta","hAllL1Eta",100,-5.0,5.0);
	hAllL1Phi = fs->make<TH1D>("hAllL1Phi","hAllL1Phi",70,-3.5,3.5);
	hAllL1Pt = fs->make<TH1D>("hAllL1Pt","hAllL1Pt",1000,0,1000.0);

	hAllL2Eta = fs->make<TH1D>("hAllL2Eta","hAllL2Eta",100,-5.0,5.0);
	hAllL2Phi = fs->make<TH1D>("hAllL2Phi","hAllL2Phi",70,-3.5,3.5);
	hAllL2Pt = fs->make<TH1D>("hAllL2Pt","hAllL2Pt",1000,0,1000.0);

	hPassL2Eta = fs->make<TH1D>("hPassL2Eta","hPassL2Eta",100,-5.0,5.0);
	hPassL2Phi = fs->make<TH1D>("hPassL2Phi","hPassL2Phi",70,-3.5,3.5);
	hPassL2Pt = fs->make<TH1D>("hPassL2Pt","hPassL2Pt",1000,0,1000.0);

	hNumGenPerEvent = fs->make<TH1D>("hNumGenPerEvent","hNumGenPerEvent",10,0.0,10.0);

	hNumL1MatchedToGen = fs->make<TH1D>("hNumL1MatchedToGen","hNumL1MatchedToGen",5,0.0,5.0);
	hNumL2MatchedToGen = fs->make<TH1D>("hNumL2MatchedToGen","hNumL2MatchedToGen",5,0.0,5.0);
	hNumL3MatchedToGen = fs->make<TH1D>("hNumL3MatchedToGen","hNumL3MatchedToGen",5,0.0,5.0);
	hNumAllL3MatchedToGen = fs->make<TH1D>("hNumAllL3MatchedToGen","hNumAllL3MatchedToGen",5,0.0,5.0);
	hNumL1MatchedToGenInEvent = fs->make<TH1D>("hNumL1MatchedToGenInEvent","hNumL1MatchedToGenInEvent",5,0.0,5.0);
	hNumL2MatchedToGenInEvent = fs->make<TH1D>("hNumL2MatchedToGenInEvent","hNumL2MatchedToGenInEvent",5,0.0,5.0);
	hNumL3MatchedToGenInEvent = fs->make<TH1D>("hNumL3MatchedToGenInEvent","hNumL3MatchedToGenInEvent",5,0.0,5.0);
	hNumAllL3MatchedToGenInEvent = fs->make<TH1D>("hNumAllL3MatchedToGenInEvent","hNumAllL3MatchedToGenInEvent",5,0.0,5.0);

	hL1Eta = fs->make<TH1D>("hL1Eta","hL1Eta",100,-5.0,5.0);
	hL1Pt = fs->make<TH1D>("hL1Pt","hL1Pt",1000,0,1000.0);
	hL2Eta = fs->make<TH1D>("hL2Eta","hL2Eta",100,-5.0,5.0);
	hL2Pt = fs->make<TH1D>("hL2Pt","hL2Pt",1000,0,1000.0);

	hL1DrGen  = fs->make<TH1D>("hL1DrGen","hL1DrGen",1500,0,1.5);
	hL2DrGen  = fs->make<TH1D>("hL2DrGen","hL2DrGen",1500,0,1.5);
	hL3DrGen  = fs->make<TH1D>("hL3DrGen","hL3DrGen",1500,0,1.5);

	hNumL2Found = fs->make<TH1D>("hNumL2Found","hNumL2Found",10,0,10.0);
	hNumL3Found = fs->make<TH1D>("hNumL3Found","hNumL3Found",10,0,10.0);
	hNumNoL3Found = fs->make<TH1D>("hNumNoL3Found","hNumNoL3Found",10,0,10.0);

	hAllDoubleL2Eta = fs->make<TH1D>("hAllDoubleL2Eta","hAllDoubleL2Eta",100,-5.0,5.0);
	hPassDoubleL2Eta = fs->make<TH1D>("hPassDoubleL2Eta","hPassDoubleL2Eta",100,-5.0,5.0);
	hAllDoubleL1Eta = fs->make<TH1D>("hAllDoubleL1Eta","hAllDoubleL1Eta",100,-5.0,5.0);
	hPassDoubleL1Eta = fs->make<TH1D>("hPassDoubleL1Eta","hPassDoubleL1Eta",100,-5.0,5.0);

	hAllDoubleL2Pt = fs->make<TH1D>("hAllDoubleL2Pt","hAllDoubleL2Pt",1000,0,1000.0);
	hPassDoubleL2Pt = fs->make<TH1D>("hPassDoubleL2Pt","hPassDoubleL2Pt",1000,0,1000.0);
	hAllDoubleL1Pt = fs->make<TH1D>("hAllDoubleL1Pt","hAllDoubleL1Pt",1000,0,1000.0);
	hPassDoubleL1Pt = fs->make<TH1D>("hPassDoubleL1Pt","hPassDoubleL1Pt",1000,0,1000.0);

}

MuonHLTAnalyser::~MuonHLTAnalyser(){}

void MuonHLTAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	std::cout << std::endl << "========== Analysing EVENT: " << iEvent.id() << "====" << std::endl;

	edm::Handle <reco::GenParticleCollection> genP;					iEvent.getByLabel("genParticles", genP);

//	Get the Reconstructed object collections:
	edm::Handle <l1extra::L1MuonParticleCollection> l1Muons;		iEvent.getByLabel(edm::InputTag("hltL1extraParticles","","BRSHLT"),l1Muons);
	edm::Handle <reco::RecoChargedCandidateCollection> l2Muons;		iEvent.getByToken(l2Token_,l2Muons);
	edm::Handle <reco::RecoChargedCandidateCollection> l3Muons;		iEvent.getByLabel(edm::InputTag("hltL3MuonCandidates","","BRSHLT"), l3Muons);

	edm::Handle<edm::View<reco::Track> > hltL3TkTracks;	iEvent.getByLabel(edm::InputTag("hltL3TkTracksFromL2OIState","","BRSHLT"), hltL3TkTracks);


//	########################### Trigger Info ###########################
	edm::Handle<edm::TriggerResults> triggerResults;		iEvent.getByLabel(edm::InputTag("TriggerResults", "", "BRSHLT"), triggerResults);
	edm::Handle<trigger::TriggerEvent> triggerSummary;		iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", "BRSHLT"), triggerSummary);

//	Get filter objects, these are the names of the paths for the Mu50 path:
	size_t L1MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag("hltL1sL1SingleMu16ORSingleMu25", "", "BRSHLT"));	//The L1 Filter
	size_t L2MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag("hltL2fL1sMu16orMu25L1f0L2Filtered16Q", "", "BRSHLT"));	//The L2 Filter
//	size_t L3MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag("hltL3fL1sMu16orMu25L1f0L2f16QL3Filtered50Q", "", "BRSHLT"));      //The L3 Filter
//	or for Mu40:
	size_t L3MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag("hltL3fL1sMu16orMu25L1f0L2f16QL3Filtered40Q", "", "BRSHLT"));      //The L3 Filter

	trigger::TriggerObjectCollection L1MuonTrigObjects;
	trigger::TriggerObjectCollection L2MuonTrigObjects;
	trigger::TriggerObjectCollection L3MuonTrigObjects;
	trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();

	if (L1MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
		//save the trigger objects corresponding to muon leg
		const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L1MuonFilterIndex);
		for (size_t j = 0; j < keysMuons.size(); j++) {
			trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
			L1MuonTrigObjects.push_back(foundObject);
		}
	}
	if (L2MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
		//save the trigger objects corresponding to muon leg
		const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L2MuonFilterIndex);
		for (size_t j = 0; j < keysMuons.size(); j++) {
			trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
			L2MuonTrigObjects.push_back(foundObject);
		}
	}
	if (L3MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
//		save the trigger objects corresponding to muon leg
		const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L3MuonFilterIndex);
		for (size_t j = 0; j < keysMuons.size(); j++) {
			trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
			L3MuonTrigObjects.push_back(foundObject);
		}
	}


	std::cout << "Number of L1s passing filter = " << L1MuonTrigObjects.size() << std::endl;
	std::cout << "Number of L2s passing filter = " << L2MuonTrigObjects.size() << std::endl;
	std::cout << "Number of L3s passing filter = " << L3MuonTrigObjects.size() << std::endl;

	std::cout << "Num L1: " << l1Muons->size() << std::endl;
	try{
		std::cout << "Num L2: " << l2Muons->size() << std::endl;
	}
	catch(...){
		std::cout << "No L2!" << std::endl;
	}
	try{
		std::cout << "Num L3: " << l3Muons->size() << std::endl;
	}
	catch(...){
		std::cout << "No L3!" << std::endl;
	}
	if (L3MuonTrigObjects.size() < L2MuonTrigObjects.size()){
		std::cout << "The Number of filter objects in L3 is less than L2 in: " << iEvent.id() << std::endl;
	}

//	################ Match the L2s (passing the L2 Filter) and L3s to Gen ################
	double NumL1MatchedToGenInEvent=0;
	double NumL2MatchedToGenInEvent=0;
	double NumL3MatchedToGenInEvent=0;
	double NumAllL3MatchedToGenInEvent=0;

	std::vector<const reco::RecoChargedCandidate*> foundL3Muons;
	std::vector<const trigger::TriggerObject*> foundL3Trig;
	std::vector<const reco::RecoChargedCandidate*> foundL3MuonsL1;

	std::vector<const reco::GenParticle*> L1FoundGens;
	std::vector<const reco::GenParticle*> L2FoundGens;
	std::vector<const reco::GenParticle*> L3FoundGens;
	std::vector<const reco::GenParticle*> L3FoundWrtL1Gens;

	int numGenPerEvent=0;

	for (unsigned int g(0); g < genP->size(); ++g){
		const reco::GenParticle* gen = &genP->at(g);
		if (fabs(gen->pdgId())==13 && gen->pt()>20 && gen->status()==1 && fabs(gen->eta())<2.4){
			++numGenPerEvent;
			hGenEta->Fill(gen->eta());
			hGenPt->Fill(gen->pt());
			hGenPhi->Fill(gen->phi());
			std::cout << "gen muon found: pt: " << gen->pt() << " eta: " << gen->eta() << " phi: " << gen->phi() << " status: " << gen->status() << std::endl;
			double NumL1MatchedToGen=0;
			double NumL2MatchedToGen=0;
			double NumL3MatchedToGen=0;
			double NumAllL3MatchedToGen=0;

			int numL2Found=0;
			int numL3Found=0;
			int numNoL3Found=0;

//			L1 Match:
			bool foundL1=false;
			for (unsigned int t(0); t < L1MuonTrigObjects.size(); ++t){
				trigger::TriggerObject* l1mu = &L1MuonTrigObjects.at(t);
				std::cout << "\tL1["<<t<<"]: deltaR(*gen,*l1mu): " << deltaR(*gen,*l1mu) << std::endl;
				hL1DrGen->Fill(deltaR(*gen,*l1mu));
				hL1Pt->Fill(l1mu->pt());
				hL1Eta->Fill(l1mu->eta());
				if (deltaR(*gen,*l1mu)<0.4 && !foundL1){
					std::cout << "\t\tL1 found: pt: " << l1mu->pt() << " eta: " << l1mu->eta() << std::endl;
					hAllL1Pt->Fill(gen->pt());
					if (gen->pt()>5){
						hAllL1Eta->Fill(gen->eta());
						hAllL1Phi->Fill(gen->phi());
						++NumL1MatchedToGen;
						++NumL1MatchedToGenInEvent;
						foundL1=true;
						L1FoundGens.push_back(gen);
					}
				}
			}

//			L2 Match:
			bool foundL2=false;
			for (unsigned int t(0); t < L2MuonTrigObjects.size(); ++t){
				trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(t);
				std::cout << "\tL2["<<t<<"] deltaR(*gen,*l2mu): " << deltaR(*gen,*l2mu) << std::endl;
				hL2DrGen->Fill(deltaR(*gen,*l2mu));
				hL2Pt->Fill(l2mu->pt());
				hL2Eta->Fill(l2mu->eta());
				if (deltaR(*gen,*l2mu)<0.05 && !foundL2){
					std::cout << "\t\tL2 found: pt: " << l2mu->pt() << " eta: " << l2mu->eta() << std::endl;
					hAllL2Pt->Fill(gen->pt());
					if (gen->pt()>5){
						hAllL2Eta->Fill(gen->eta());
						hAllL2Phi->Fill(gen->phi());
						foundL2=true;
						++NumL2MatchedToGen;
						++NumL2MatchedToGenInEvent;
						L2FoundGens.push_back(gen);
					}
				}
			}
//			L3 Match:
			if (foundL2){
				++numL2Found;
				bool foundL3=false;

				for (unsigned int t(0); t < L3MuonTrigObjects.size(); ++t){
					trigger::TriggerObject* trigMuObj = &L3MuonTrigObjects.at(t);
					hL3DrGen->Fill(deltaR(*gen,*trigMuObj));
					std::cout << "\tL3["<<t<<"]: deltaR(*gen,*l3mu): " << deltaR(*gen,*trigMuObj) << std::endl;
					if (deltaR(*gen,*trigMuObj)<0.01 && !foundL3){
						if (std::find(foundL3Trig.begin(), foundL3Trig.end(), trigMuObj)!=foundL3Trig.end()) std::cout << "THIS L3 WAS ALREADY FOUND!" << std::endl;
						foundL3Trig.push_back(trigMuObj);
						std::cout << "\t\tL3 found: pt: " << trigMuObj->pt() << " eta: " << trigMuObj->eta() << std::endl;
						hPassL2Pt->Fill(gen->pt());
						if (gen->pt()>5){
								hPassL2Eta->Fill(gen->eta());
								hPassL2Phi->Fill(gen->phi());
								++NumL3MatchedToGen;
								++NumL3MatchedToGenInEvent;
								L3FoundGens.push_back(gen);
								foundL3=true;
						}
					}
				}

			} //FoundL2

			hNumL1MatchedToGen->Fill(NumL1MatchedToGen);
			hNumL2MatchedToGen->Fill(NumL2MatchedToGen);
			hNumL3MatchedToGen->Fill(NumL3MatchedToGen);
			hNumAllL3MatchedToGen->Fill(NumAllL3MatchedToGen);

			hNumL2Found->Fill(numL2Found);
			hNumL3Found->Fill(numL3Found);
			hNumNoL3Found->Fill(numNoL3Found);
			std::cout << "numL2Found: " << numL2Found
					<< " numL3Found: " << numL3Found
					<< " numNoL3Found: " << numNoL3Found
					<< std::endl;

		} //genMuon
	} //genParticle col

	std::cout << "Filling Double Muon Efficiency: NumGen = " << numGenPerEvent
			<< "\n\tL2FoundGens.size(): " << L2FoundGens.size()
			<< " L3FoundGens.size(): " << L3FoundGens.size()
			<< "\n\tL1FoundGens.size(): " << L1FoundGens.size()
			<< " L3FoundWrtL1Gens.size(): " << L3FoundWrtL1Gens.size()
			<< std::endl;


	if (L2FoundGens.size()==2){
		std::cout << "Yes there are 2 L2FoundGens: " << L2FoundGens.size() << std::endl;
		const reco::GenParticle* gen = L2FoundGens.at(0);
		hAllDoubleL2Eta->Fill(gen->eta());
		hAllDoubleL2Pt->Fill(gen->pt());
		if (L3FoundGens.size()==2){
			hPassDoubleL2Eta->Fill(gen->eta());
			hPassDoubleL2Pt->Fill(gen->pt());
		}
		const reco::GenParticle* gen2 = L2FoundGens.at(1);
		hAllDoubleL2Eta->Fill(gen2->eta());
		hAllDoubleL2Pt->Fill(gen2->pt());
		if (L3FoundGens.size()==2){
			hPassDoubleL2Eta->Fill(gen2->eta());
			hPassDoubleL2Pt->Fill(gen2->pt());
		}

	}
	if (L1FoundGens.size()==2){
		std::cout << "Yes there are 2 L1FoundGens: " << L1FoundGens.size() << std::endl;
		const reco::GenParticle* gen = L1FoundGens.at(0);
		hAllDoubleL1Eta->Fill(gen->eta());
		hAllDoubleL1Pt->Fill(gen->pt());
		if (L3FoundWrtL1Gens.size()==2){
			hPassDoubleL1Eta->Fill(gen->eta());
			hPassDoubleL1Pt->Fill(gen->pt());
		}
		const reco::GenParticle* gen2 = L1FoundGens.at(1);
		hAllDoubleL1Eta->Fill(gen2->eta());
		hAllDoubleL1Pt->Fill(gen2->pt());
		if (L3FoundWrtL1Gens.size()==2){
			hPassDoubleL1Eta->Fill(gen2->eta());
			hPassDoubleL1Pt->Fill(gen2->pt());
		}
	}
	hNumGenPerEvent->Fill(numGenPerEvent);


	hNumL1MatchedToGenInEvent->Fill(NumL1MatchedToGenInEvent);
	hNumL2MatchedToGenInEvent->Fill(NumL2MatchedToGenInEvent);
	hNumL3MatchedToGenInEvent->Fill(NumL3MatchedToGenInEvent);
	hNumAllL3MatchedToGenInEvent->Fill(NumAllL3MatchedToGenInEvent);
	std::cout << "NumL1MatchedToGenInEvent: " << NumL1MatchedToGenInEvent
							<< " NumL2MatchedToGenInEvent: " << NumL2MatchedToGenInEvent
							<< " NumL3MatchedToGenInEvent: " << NumL3MatchedToGenInEvent
							<< " NumAllL3MatchedToGenInEvent: " << NumAllL3MatchedToGenInEvent
							<< std::endl;

	if (NumL3MatchedToGenInEvent<NumL2MatchedToGenInEvent) std::cout << "L3 Not Found!: " << iEvent.id().luminosityBlock() << ":" << iEvent.id().event() << std::endl;


//	##################Plot seed and tracks information#################

//	Seeds and Tracker Tracks:
	try{
		try{
			edm::Handle<L3MuonTrajectorySeedCollection> hltL3TrajSeedOIState;
			iEvent.getByLabel(edm::InputTag("hltL3TrajSeedOIState","","BRSHLT"), hltL3TrajSeedOIState);
			histhltL3TrajSeedOIState->Fill(hltL3TrajSeedOIState->size());
			std::cout << "# of hltL3TrajSeedOIState: " << hltL3TrajSeedOIState->size() << std::endl;

			edm::Handle<edm::View<reco::Track> > hltL3TkTracksFromL2OIState;
			iEvent.getByLabel(edm::InputTag("hltL3TkTracksFromL2OIState","","BRSHLT"), hltL3TkTracksFromL2OIState);
			histhltL3TkTracksFromL2OIState->Fill(hltL3TkTracksFromL2OIState->size());
			std::cout << "# of hltL3TkTracksFromL2OIState = " << hltL3TkTracksFromL2OIState->size() << std::endl;
		}
		catch(...){
			std::cout << "Could not get OIState info!" << std::endl;
		}

		try{
			edm::Handle<L3MuonTrajectorySeedCollection> hltL3TrajSeedOIHit;
			iEvent.getByLabel(edm::InputTag("hltL3TrajSeedOIHit","","BRSHLT"), hltL3TrajSeedOIHit);
			histhltL3TrajSeedOIHit->Fill(hltL3TrajSeedOIHit->size());
			std::cout << "# of hltL3TrajSeedOIHit: " << hltL3TrajSeedOIHit->size() << std::endl;

			edm::Handle<edm::View<reco::Track> > hltL3TkTracksFromL2OIHit;
			iEvent.getByLabel(edm::InputTag("hltL3TkTracksFromL2OIHit","","BRSHLT"), hltL3TkTracksFromL2OIHit);
			histhltL3TkTracksFromL2OIHit->Fill(hltL3TkTracksFromL2OIHit->size());
			std::cout << "# of hltL3TkTracksFromL2OIHit = " << hltL3TkTracksFromL2OIHit->size() << std::endl;
		}
		catch(...){
			std::cout << "Could not get OIHit info!" << std::endl;
		}

		try{
			edm::Handle<L3MuonTrajectorySeedCollection> hltL3TrajSeedIOHit;
			iEvent.getByLabel(edm::InputTag("hltL3TrajSeedIOHit","","BRSHLT"), hltL3TrajSeedIOHit);
			histhltL3TrajSeedIOHit->Fill(hltL3TrajSeedIOHit->size());
			std::cout << "# of hltL3TrajSeedIOHit: " << hltL3TrajSeedIOHit->size() << std::endl;

			edm::Handle<edm::View<reco::Track> > hltL3TkTracksFromL2IOHit;
			iEvent.getByLabel(edm::InputTag("hltL3TkTracksFromL2IOHit","","BRSHLT"), hltL3TkTracksFromL2IOHit);
			histhltL3TkTracksFromL2IOHit->Fill(hltL3TkTracksFromL2IOHit->size());
			std::cout << "# of hltL3TkTracksFromL2IOHit = " << hltL3TkTracksFromL2IOHit->size() << std::endl;
		}
		catch(...){
			std::cout << "Could not get IOHit info!" << std::endl;
		}
	}
	catch(...){
		std::cout << "Could not get seeds/tracker tracks information for L3 Cascade!!" << std::endl;
	}

	std::cout << "============================================================" << std::endl
			<< "============================================================" << std::endl;

	std::cout << std::endl<< std::endl<< std::endl<< std::endl<< std::endl;

}

void MuonHLTAnalyser::beginJob(){}

void MuonHLTAnalyser::endJob(){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonHLTAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonHLTAnalyser);
