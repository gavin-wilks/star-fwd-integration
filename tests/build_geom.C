//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable


void build_geom( TString geomtag = "dev2021") {
   gROOT->LoadMacro("bfc.C");
   bfc(0, "fzin agml sdt20181215", "" );

   gSystem->Load("libStarClassLibrary.so");
   gSystem->Load("libStEvent.so" );

   // Force build of the geometry
   TFile *geom = TFile::Open( "fGeom_"+geomtag+".root" );

   if ( 0 == geom ) {
      AgModule::SetStacker( new StarTGeoStacker() );
      AgPosition::SetDebug(2);
      StarGeometry::Construct(geomtag.Data());

      // Genfit requires the geometry is cached in a ROOT file
      gGeoManager->Export( "fGeom_"+geomtag+".root" );
      cout << "Writing output to geometry file [" << "fGeom_" << geomtag.Data() << ".root" << "]" << endl;
   }
   else {
      cout << "WARNING:  Geometry file [" << "fGeom_" << geomtag.Data() << ".root" << "] already exists." << endl;
      cout << "Existting without doing anything!" << endl;
      delete geom;
   }

}
