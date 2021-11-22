//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

TFile *output = 0;

void slow_track(  int n = 100,
                  const char *inFile = "out/sim_2.fzd",
                  std::string configFile = "slow_track.xml",
                  const char *geom = "dev2021",
                  std::string qaoutname = "") {
    TString _geom = geom;

    //bool SiIneff = false;    

    TString _chain;

    // NOTE "event" does not work in CMAKE StRoot wo network, it includes detDb - root problem. Swap to StEvent instead
    _chain = Form("fzin %s StEvent evout geantout ReverseField agml usexgeom bigbig", _geom.Data());

    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    // StarMagField::setConstBz(true);

    gSystem->Load("libMathMore.so");
    gSystem->Load("libXMLIO.so"); // needed by FwdTrackerConfig
    gSystem->Load("libStarGeneratorUtil.so"); // needed for StarRandom
    gSystem->Load("libStFstSimMaker.so");
    gSystem->Load("libStFttSimMaker.so");

    gSystem->Load("libgenfit2.so"); // needed for GenFit
    gSystem->Load("libKiTrack.so"); // needed for KiTrack
    gSystem->Load("libStEventUtilities.so");
    gSystem->Load("libStFwdTrackMaker.so");

    StFttSlowSimMaker *fttSlowSim = new StFttSlowSimMaker();
    cout << "Adding StFttSlowSimMaker to chain" << endl;
    chain->AddMaker(fttSlowSim);

    // Create fast simulator and add after event maker
    StFstSlowSimMaker *fstSlowSim = new StFstSlowSimMaker();

    fstSlowSim->SetInEfficiency(true, 0.05, 0.10); // inefficiency of Si 
    fstSlowSim->SetCrossTalk(true);    

    //!!!!!   for mc hits   !!!!!//
    fstSlowSim->SetFillHist( true );
    
    fstSlowSim->SetQAFileName(qaoutname.c_str());

    cout << "Adding StFstSlowSimMaker to chain" << endl;
    chain->AddMaker(fstSlowSim);

    /*StFwdTrackMaker *gmk = new StFwdTrackMaker();
    // config file set here overides chain opt
    gmk->SetConfigFile( configFile );
    gmk->SetGenerateTree( true );
    gmk->SetGenerateHistograms( true ); //!!!!!   for mc hits   !!!!!// 
    chain->AddMaker(gmk);
*/
    chain->Init();

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {

        chain->Clear();
        if (kStOK != chain->Make())
            break;
    }
     
}
