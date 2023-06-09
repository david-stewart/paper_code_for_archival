#include <TSystem>
#include <ctime>
#include <string>

class StMaker;
class StChain;
class StPicoDstMaker;

StChain *chain;
void myExec(
    const Int_t   nEv_input=-1,
    const Char_t *picoDstList="short_pico_union.list",
    const Char_t *out_name   ="test",
    const char*   bad_tow_list ="bad_tower.list",
    const char*   bad_run_list ="bad_run.list",
    const int     debug_level=0
){
    //------------------------------------------------------
    // Load required libraries
    //------------------------------------------------------
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    gSystem->Load("StMuDSTMaker");  
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StChain");
    gSystem->Load("StDbBroker");
    gSystem->Load("St_db_Maker.so");
    gSystem->Load("StEmcUtil");
    gSystem->Load("$FASTJET/lib/libfastjet");
    gSystem->Load("$FASTJET/lib/libfastjettools");

    // load local (and common) libraries from ./StRoot
    gSystem->Load("TreeObj");
	gSystem->Load("JetBranchMaker");
    gSystem->Load("MemTimeProgression");
    gSystem->Load("IntList");
    gSystem->Load("TrackTowerMatcher");
    gSystem->Load("SL20a_picoTrackTowerMaker");
    gSystem->Load("myExec");

    //-----------------------------------------
    //  output file
    //-----------------------------------------
    ofstream log;
    log.open(Form("%s.log",out_name));
    if (!log.is_open()) {
        cout << "Fatal error: could not open output file " << Form("%s.root",out_name) << endl;
        return kFatal;
		}

    //-----------------------------------------
    //  get run parameters of the file list
    //-----------------------------------------
    log << "Reading input picoDst list: " << picoDstList << endl;

    const char* database = "mysql://db04.star.bnl.gov:3414/RunLog?timeout=60";
    const char* user = "djs232";
    const char* pass = "";
    TMySQLServer* mysql = TMySQLServer::Connect(database,user,pass);
    
    ifstream file;
    file.open(picoDstList);
    if (!file.is_open()){
        cout << "Fatal error: Could not open " << picoDstList << endl;
        log  << "Fatal error: Could not open " << picoDstList << endl;
        return kFatal;
    }
    string str;
    ifstream file;
    file.close();
    mysql->Close();
    log << "Finished run ids, times, and durations. Closing mysql database." << endl << endl;

    //-----------------------------------------
    // process the input file for xrootd files
    //-----------------------------------------
    // process in the input file to get rid of the run numbers (if present)
    // for when the input is like: 
    // /tmp/djs232/D96D05E45C742BAAF8F88D9560575DEE_0/INPUTFILES/st_ssdmb_16124017_raw_2500003.picoDst.root 8586
    // but the picoDst wants the line without the 8586 appended
    file.open(picoDstList);
    ofstream fout;
    fout.open("linted_picos.list");
    while(getline(file,str)){
        TString tstr = str;
        cout << "Line before: " << tstr << endl;
        TPRegexp("(\\S+) (\\d+)").Substitute(tstr,("$1"));
        fout << tstr << endl;
        cout << "Line after " << tstr << endl;
    };
    fout.close();

    //----------------------------------------
    // Get the bemc maker
    //----------------------------------------
    chain = new StChain();
    log << "Opening BEMC database." << endl;


    //-----------------------------------
    // Make chain and run main macro
    //-----------------------------------
    StPicoDstMaker *picoMaker =0X0; 
    StPicoDstMaker::PicoIoMode IoMode = 2;
    cout << "picoMaker0: " << picoMaker << endl;
    cout << "picoDstList: " << picoDstList << endl;
    picoMaker = new StPicoDstMaker(IoMode,"linted_picos.list","picoDst");
    picoMaker->SetStatus("*",1);

    cout << "picoMaker1: " << picoMaker << endl;
    /* picoMaker->SetBranchStatus("*",1); */


    St_db_Maker* starDb = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb");
    starDb->SetFlavor("ofl");
    starDb->SetAttr("blacklist", "eemc");
    starDb->SetAttr("blacklist", "svt");
    starDb->SetAttr("blacklist", "tpc");
    starDb->SetAttr("blacklist", "tof");
    starDb->SetAttr("blacklist", "pxl");
    starDb->SetAttr("blacklist", "ist");
    starDb->SetAttr("blacklist", "bprs");
    starDb->SetAttr("blacklist", "smd");
    starDb->SetAttr("blacklist", "bsmd");

    myExec *myrunQA = new myExec(
        "pf",
        out_name,
        log,
        picoMaker,
        bad_run_list,
        bad_tow_list,
        starDb,
        debug_level
    );

    log << "chain->Init()" <<endl;
    chain->Init();

    int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    log  << "- Total events: = " << total << endl;
    int nEvents = nEv_input;
    if(nEvents>total||nEvents<0) nEvents = total;
    cout<<"Number of Events = "<<nEvents<<endl;
    log <<"- Number of events to read = "<<nEvents<<endl;

    int usage0 = 0;
    for (Int_t i=0; i<nEvents; i++){   
        if(i%5000==0) {
            cout << "- Finished " << i << " events." << endl;
            log  << "- Finished " << i << " events." << endl;

        }
        chain->Clear();
        int iret = chain->Make(i);
        /* log << " i: " << i << endl; */
        if (iret) { cout << "Bad return code!" << iret << endl; break;}
        if (iret) { log  << "Bad return code!" << iret << endl; break;}
        total++;
    } 

    current_time = time(NULL);
    loc = localtime(&current_time);

	log << endl;	
	log << "----------------------" << endl;
	log << "|  Script complete!  |"<< endl;
	log << "----------------------" << endl;
	chain->Finish();

	delete chain;
}
