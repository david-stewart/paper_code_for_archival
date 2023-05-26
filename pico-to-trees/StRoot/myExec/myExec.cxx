//   New objective: kj
// Includes code from (1) Joel Mazer from https://github.com/joelmazer/star-jetframework in July 2018
//                    (2) Tong Liu from ARCAS, also July 2018

#include "myExec.h"
#include <iostream>
#include "TChain.h"
#include <fstream>
#include "TNtuple.h"
#include "TF1.h"
#include "fastjet/config.h"
#include <algorithm>
#include "TrackTowerMatcher/TrackTowerMatcher.h"
#include "GetMemValue.h"

#define  MIN_TRACK_PT 0.2
#define  MIN_TOWER_ET 0.2
#define  MaxNTracks 350
#define  MaxNTowers 700 // above 4 GeV
/* #define  MaxNJets   40 // above 4 GeV */

ClassImp(myExec)
//-----------------------------------------------------------------------------
myExec::myExec(
        const char*     _name,
        const char*     _outName,
        ofstream&       _log,
        StPicoDstMaker* _picoMaker,
        const char*     _bad_run_list,
        const char*     _bad_tower_list,
        St_db_Maker*    _starDb,
        int             _debug_level
    )
  : StMaker         (_name), 
    log             (_log),
    picoDstMaker    {_picoMaker},
	bad_run_list    {_bad_run_list,log},
	bad_tower_list  {_bad_tower_list,log},

    starDb            {_starDb},
    debug_level       {_debug_level},
	stats             {500},

    bemcGeom          { StEmcGeom::getEmcGeom("bemc") },
    fout_root         { new TFile(Form("%s.root",_outName),"RECREATE") },
    ftree             { new TTree("events","tree many triggers") },
    picoTT            { _picoMaker, bad_tower_list, log, starDb,
                        _debug_level, MaxNTracks, MaxNTowers, false},
    clones_mupicoTrack          { "mupicoTrack", MaxNTracks },
    clones_mupicoTower          { "mupicoTower", MaxNTowers },
	/* C_jets_maker      { ftree, "Cjet", true, "anti-kt", MaxNJets }, */
	/* F_jets_maker      { ftree, "Fjet", true, "anti-kt", MaxNJets }, */ 
    // add triggers: 55, 500018, 500202, 32, 500904, 57
    ftriggers         { ftree, 
        {500004, 500008, 500018, 500205, 500206, 500215, 500001,
            55, 500202, 32, 500904, 57
        }}
{
    ftree->Branch("event",   &fevent);
    ftree->Branch("track",   &clones_mupicoTrack);
    ftree->Branch("tower",   &clones_mupicoTower);

    ftree->Branch("BbcEastRate",   &BbcEastRate);
    ftree->Branch("BbcWestRate",   &BbcWestRate);
    ftree->Branch("nGlobal",   &nGlobal);

	/* ftree->Branch("EastBBC", east, "s[24]/s"); */
	/* ftree->Branch("WestBBC", west, "s[24]/s"); */

	/* ftree->Branch("ZdcSmdEastHorizontal", ZdcSmdEastHorizontal, "s[8]/s"); */
	/* ftree->Branch("ZdcSmdEastVertical", ZdcSmdEastVertical, "s[8]/s"); */
	/* ftree->Branch("ZdcSmdWestHorizontal", ZdcSmdWestHorizontal, "s[8]/s"); */
	/* ftree->Branch("ZdcSmdWestVertical", ZdcSmdWestVertical, "s[8]/s"); */
};

//----------------------------------------------------------------------------- 
myExec::~myExec()
{

};

//----------------------------------------------------------------------------- 
Int_t myExec::Init() {
    return kStOK;
};

//----------------------------------------------------------------------------- 
Int_t myExec::Finish() {
    log << " Finished with myExec object " << endl;
    cout << " Finished with myExec object " << endl;
    cout << "----------------------------------------------------"<< endl;
    stats.set_get_stats();
    log  << " Final stats: " << stats.stats << endl;
    cout << " Max memory used: " << stats.max_mem/1000. << " MB " << endl;
    log  << " Max memory used: " << stats.max_mem/1000. << " MB " << endl;
    cout << "----------------------------------------------------"<< endl;
	/* ftree->Write(); */
    fout_root->Write();
    fout_root->Close();
    return kStOK;
};
    
//----------------------------------------------------------------------------- 
void myExec::Clear(Option_t *opt) { };
//----------------------------------------------------------------------------- 
Int_t myExec::Make() {
	if (stats.call()) log << stats.stats << endl;
    if(!picoDstMaker) {
        log << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    StPicoDst* picoDst = picoDstMaker->picoDst();
    if(!picoDst) {
        log << " No PicoDst! Skip! " << endl;
        return kStWarn;
    }

    StPicoEvent* picoEvent = picoDst->event();

    for (auto trig : ftriggers.trigger_ids) {
        ftriggers.set_trigger(trig, picoEvent->isTrigger(trig));
    }
    if (! ftriggers.has_any_triggers() ) return kStOk;

    fevent.refMult = picoEvent->refMult();
    fevent.runId = picoEvent->runId() ;
    fevent.ZDCx = picoEvent->ZDCx();
    fevent.BBCx = picoEvent->BBCx();
    fevent.eventId = picoEvent->eventId();

    nGlobal = picoEvent->numberOfGlobalTracks();
    BbcEastRate = picoEvent->bbcEastRate();
    BbcWestRate = picoEvent->bbcWestRate();

    if (debug_level > 3) cout << " debug:  runId " << fevent.runId << endl;

	if (bad_run_list(fevent.runId)) return kStOk;

    // cut on vz
    fevent.vx = picoEvent->primaryVertex().x();
    fevent.vy = picoEvent->primaryVertex().y();
    fevent.vz = picoEvent->primaryVertex().z();
    fevent.vzVpd = picoEvent->vzVpd();
    if (TMath::Abs(fevent.vz) > 100.) return kStOk;
    /* if (TMath::Abs(fevent.vz-fevent.vzVpd) > 10.) return kStOk; */

    fevent.ranking = picoEvent->ranking();
    if (fevent.ranking < 0) return kStOk;

    clones_mupicoTrack.Clear();
    clones_mupicoTower.Clear();
	/* C_jets_maker.clear(); */
	/* F_jets_maker.clear(); */

    picoTT.make();

    // Fill in tracks (and charged jets) from picoDst
    int i{0};
    /* cout << " here " << endl; */
    while (picoTT.next_track()){
         mupicoTrack* track = (mupicoTrack*) clones_mupicoTrack.ConstructedAt(i);
         *track = picoTT.track();
            /* cout << Form(" nHitsFit:  %i     nHitsPoss:  %i   pT:  %5.2f", */
                /* track->nHitsFit, track->nHitsPoss, track->pt) << endl; */
         /* if (track->pass_cuts) */
             /* C_jets_maker.fill( {{track->pt,track->eta,track->phi},i} ); */
         ++i;
    }
    /* C_jets_maker.generate(); */
    
    // Fill in towers from muDst (and into jets)
    /* F_jets_maker.fill(C_jets_maker.particles); */
    i=0;
    while (picoTT.next_tower()){
         mupicoTower* tower = (mupicoTower*) clones_mupicoTower.ConstructedAt(i);
         *tower = picoTT.tower();

         // fill in jet ; towers have been already cut for bad_tower_list
         /* float Et { tower->Et - tower->Et_hadroncorr }; */
         /* if (Et >= 0.2) F_jets_maker.fill( {{Et,tower->eta,tower->phi},-i-1} ); */
         ++i;
    }
    /* F_jets_maker.generate(); */

	
    // fill in BBC event
    for (int i = 0; i < 24; ++i) {
		fevent.EastBBC[i] = picoEvent->bbcAdcEast(i);
		fevent.WestBBC[i] = picoEvent->bbcAdcWest(i);
	}

    fevent.BBC_Ein = 0.;
    fevent.BBC_Win = 0.;
    for (int i = 0; i < 16; ++i) {
        fevent.BBC_Ein += fevent.EastBBC[i];
        fevent.BBC_Win += fevent.WestBBC[i];
    }
    fevent.BBC_Wout = 0.;
    fevent.BBC_Eout = 0.;
    for (int i = 16; i<24; ++i) {
        fevent.BBC_Eout += fevent.EastBBC[i];
        fevent.BBC_Wout += fevent.WestBBC[i];
    }

    fevent.ZdcSumAdcEast = picoEvent->ZdcSumAdcEast();
    fevent.ZdcSumAdcWest = picoEvent->ZdcSumAdcWest();

    for (int i=0; i<8; ++i) {
        fevent.ZdcSmdWestVertical[i]   = picoEvent->ZdcSmdWestVertical(i);
        fevent.ZdcSmdWestHorizontal[i] = picoEvent->ZdcSmdWestHorizontal(i);
        fevent.ZdcSmdEastVertical[i]   = picoEvent->ZdcSmdEastVertical(i);
        fevent.ZdcSmdEastHorizontal[i] = picoEvent->ZdcSmdEastHorizontal(i);
    }

    if (debug_level == 3) cout << "debug: filling" << endl;
    ftree->Fill();
    return kStOk;
};
