// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
// root4star starsim.C
// This macro is capable of loading the misalignment for the FST geometry using dev2022m geomtag
// To ensure the misalignment tables load, the SDT needs to be set later than timestamp on the tables
// 	located under ./StarDb/Geometry/fst/ 
class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary  = 0;

class StarKinematics;
StarKinematics *kinematics  = 0;

int     _npart = 10;    // floor number of tracks per event
TString _part  = "mu+"; // particle to simulate
float   _ptmn  = 0.200; // min pT to simulate [GeV]
float   _ptmx  = 5.000; // max pT to simulate [GeV]
float   _etamn = 2.5;   // min eta to simulate
float   _etamx = 4.0;   // max eta to simulate

TString _geometry = "dev2022m";
TString DBV; 
TString SDT = "sdt20220421";

//______________________________________________________________________________________

void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag;
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}

//______________________________________________________________________________________

void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}

//______________________________________________________________________________________

void trig( int n=1 )
{


  for ( int i=0; i<n; i++ ) {

    // Clear the chain from the previous event
    chain->Clear();

    _primary->SetVertex( 0.0, 0.0, 0.0 );
    _primary->SetSigma(0.1, 0.1, 30.0);

    kinematics->Kine( _npart, _part, _ptmn, _ptmx, _etamn, _etamx );

    // Generate the event
    chain->Make();

  }
}
//______________________________________________________________________________________

void Kinematics()
{
  gSystem->Load( "libKinematics.so");
  kinematics = new StarKinematics();
    
  _primary->AddGenerator(kinematics);
}
//______________________________________________________________________________________

void starsim ( int rngSeed=0,      
               int nevents=-1,    
               const char* outfile = "out" )
{ 
 
  //
  //________________________________________________________
  //
  // Setup the big full chain
  //
  //________________________________________________________
  //
  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = _geometry; simple += " ";
    simple += SDT; simple += " ";
    simple += DBV; simple += " ";
    simple += " geant gstar usexgeom misalign agml ";
    bfc(0, simple );
  }
  //
  //________________________________________________________
  //
  // Load in supporting libraries
  //________________________________________________________
  //
  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so"  );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so"  );
  gSystem->Load( "StarGeneratorDecay.so" ); 
  gSystem->Load( "libMathMore.so"        );
  // Force xgeometry.so from local repository to load   
  gSystem->Load( "xgeometry.so"          );


  //________________________________________________________
  //
  // Setup RNG seed and map all ROOT TRandom here
  //________________________________________________________
  // 
 
  StarRandom::seed( rngSeed ); // but will reset based on run number and event number
  StarRandom::capture();
  
  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "kinematics.starsim.root");
    chain -> AddBefore( "geant", _primary );
  }

  Kinematics();

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();

  //
  // Setup geometry and set starsim to use agusread for input
  //
  command("gkine -4 0");

  TString fzdname = Form("gfile o %s_%i.fzd",outfile,rngSeed);
  TString rooname = Form("%s_%i.genevents.root",outfile,rngSeed);
  command( fzdname );
  _primary -> SetFileName( rooname );

  //
  // Trigger on nevents
  //
  trig( nevents );

  command("call agexit");  // Make sure that STARSIM exits properly

}

//______________________________________________________________________________________

void starsim ( const char* outfile,
               const int nevents,
               const int np,
               const char* part,
               const float ptmn,
               const float ptmx,
	       const float etamn,
               const float etamx,
               const int seed,
               const char* geom,
               const char* dbv = 0)
{
  _npart = np;
  _part  = part;
  _ptmn  = ptmn;
  _ptmx  = ptmx;
  _etamn = etamn;
  _etamx = etamx;
  _geometry = geom;
  if ( dbv ) DBV = dbv;
  starsim( seed, nevents, outfile );
}
//______________________________________________________________________________________

