#ifndef myExec_h
#define myExec_h

#include "TTree.h"
/* #include "St_db_Maker/St_db_Maker.h" */
/* #include "StEmcUtil/database/StBemcTables.h" */
/* #include "StRoot/StPicoDstMaker/StPicoDst.h" */
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include <iostream>
#include "StThreeVectorF.hh"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"

/* #include "TreeEvent/TreeEvent.h" */
#include "TreeObj/TreeObj.h"
#include "JetBranchMaker/JetBranchMaker.h"
#include "MemTimeProgression/MemTimeProgression.h"
#include "IntList/IntList.h"
#include "SL20a_picoTrackTowerMaker/SL20a_picoTrackTowerMaker.h"

#include <map>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class TString;

class myExec : public StMaker {

    private:
    // initialized from constructor
    ofstream&       log;
    StPicoDstMaker* picoDstMaker;
	IntList         bad_run_list;
    IntList         bad_tower_list;
    St_db_Maker*    starDb;
    int             debug_level;
	MemTimeProgression stats;
    /* bool is_trans(double, double); */

    // other
    /* vector<int>     bad_run_list; */

    vector<int>     triggers;
    StBemcTables*   bemc;
    StEmcGeom*      bemcGeom;

    vector<int>     v_triggers;
    TFile*          fout_root;
    TTree*          ftree;


    SL20a_picoTrackTowerMaker picoTT;
    TClonesArray    clones_mupicoTrack;
    TClonesArray    clones_mupicoTower;
	/* JetBranchMaker  C_jets_maker; // charged jets maker */
	/* JetBranchMaker  F_jets_maker; // full jets maker */

    mupicoEventHeader fevent;
    myTriggers        ftriggers;

    float _BBCx;
    short nGlobal;
    float BbcEastRate;
    float BbcWestRate;

    /* UShort_t east[24]; */
    /* UShort_t west[24]; */

    /* UShort_t ZdcSmdEastHorizontal[8]; */
    /* UShort_t ZdcSmdEastVertical[8]; */
    /* UShort_t ZdcSmdWestHorizontal[8]; */
    /* UShort_t ZdcSmdWestVertical[8]; */

    public:
    myExec(
        const char*     _name,
        const char*     _outName,
        ofstream&       _log,
        StPicoDstMaker* _picoMaker,
        const char*     _bad_run_list,
        const char*     _bad_tower_list,
        St_db_Maker*    _starDb,
        int             _debug_level = 0
    );
    virtual ~myExec();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    ClassDef(myExec, 1)
};

#endif
