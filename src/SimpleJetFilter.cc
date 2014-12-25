// -*- C++ -*-
//
// Package:    Utilities
// Class:      SimpleJetFilter
// 
/**\class SimpleJetFilter SimpleJetFilter.cc DPGAnalysis/Skims/src/SimpleJetFilter.cc

 Description: 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea Venturi
//         Created:  Tue Oct 21 20:55:22 CEST 2008
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/JetReco/interface/PFJet.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

//
// class declaration
//

class SimpleJetFilter : public edm::EDFilter 
{
public:
  explicit SimpleJetFilter(const edm::ParameterSet&);
  ~SimpleJetFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
      
      // ----------member data ---------------------------

  edm::InputTag m_jetCollection;
  const double m_ptcut;
  const double m_etamaxcut;
  const double m_njetmin;
  PFJetIDSelectionFunctor m_jetIDfunc;

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
SimpleJetFilter::SimpleJetFilter(const edm::ParameterSet& iConfig):
  m_jetCollection(iConfig.getParameter<edm::InputTag>("jetCollection")),
  m_ptcut(iConfig.getParameter<double>("ptCut")),
  m_etamaxcut(iConfig.getParameter<double>("maxRapidityCut")),
  m_njetmin(iConfig.getParameter<unsigned int>("nJetMin")),
  m_jetIDfunc(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::LOOSE)
{
   //now do what ever initialization is needed


}

SimpleJetFilter::~SimpleJetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SimpleJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  bool selected = false;
  
  Handle<reco::PFJetCollection> jetcoll;
  iEvent.getByLabel(m_jetCollection,jetcoll);
  
  unsigned int goodjets = 0;

  for(unsigned int ijet=0;ijet<jetcoll->size();++ijet) {
    
    const reco::PFJetRef jet(jetcoll,ijet);

    LogDebug("JetUnderTest") << "Jet with eta = " << jet->eta() << " and pt = " << jet->pt() << " under test";

    if( !(std::abs(jet->eta()) < m_etamaxcut && jet->pt() > m_ptcut )) continue;

    LogDebug("JetUnderTest") << "kincut passed";

    pat::strbitset ret = m_jetIDfunc.getBitTemplate();
    ret.set(false);
    bool goodjet = m_jetIDfunc((*jetcoll)[ijet],ret);
    if(goodjet) { 
      ++goodjets;
      LogDebug("JetUnderTest") << "JetID passed";
    }
    if(goodjets >= m_njetmin) return true;
    
  }  
  return selected;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SimpleJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleJetFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleJetFilter);
