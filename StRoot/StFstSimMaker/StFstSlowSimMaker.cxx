#include "StFstSlowSimMaker.h"

#include "St_base/StMessMgr.h"

#include "StEvent/StEvent.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarGenerator/UTIL/StarRandom.h"

#include "TCanvas.h"
#include "TCernLib.h"
#include "TH2F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector2.h"
#include "TVector3.h"

#include <array>
#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <string>


// lets not polute the global scope
namespace FstGlobalSlow{
	// converts the position error to a cov mat
    StMatrixF Hack1to6(const StHit *stHit);
	// used to convert between uniform and 1 sigma error
    constexpr float SQRT12 = sqrt(12.0);

    // Disk segmentation
    //
    float RMIN[] = {0.95 * 4.3, 0.95 * 4.3, 0.95 * 4.3, 0.95 * 5.0, 0.95 * 5.0, 0.95 * 5.0};
    float RMAX[] = {1.05 * 15.0, 1.05 * 25.0, 1.05 * 25.0, 1.05 * 28.0, 1.05 * 28.0, 1.05 * 28.0};

    //NEXT IS only for disk ARRAY 456 with the radius from 5 to 28.
    float RSegment[] = {5., 7.875, 10.75, 13.625, 16.5, 19.375, 22.25, 25.125, 28.};

    // controls some extra output
    const bool verbose = false;

	// key type for lookup tables on disk r-phi strip
	typedef std::tuple<int, int, int> FstKeyTriple;

}

StFstSlowSimMaker::StFstSlowSimMaker(const Char_t *name)
	: StMaker{name},
    mNumR{8},
    mNumPHI{128},
    mNumSEC{12},
    mRaster{0},
    mInEff{0},
    mInnerSiInEff{0},
    mOuterSiInEff{0},
    mCrossTalk{0},
    mHist{0},
    mQAFileName(0),
    hMCHit{0},
    hMCPhiZOut{0},
    hMCPhiZIn{0},
    hRCHit{0},
    hRCPhiZOut{0},
    hRCPhiZIn{0},
    hTrutHitYXDisk(0),
    hTrutHitRDisk(0),
    hTrutHitRShower{0},
    hTrutHitPhiDisk(0),
    hTrutHitPhiZ(0),
    hRecoHitYXDisk(0),
    hRecoHitRDisk(0),
    hRecoHitPhiDisk(0),
    hRecoHitPhiZ(0),
    hGlobalDRDisk(0),
    hGlobalZ(0),
    h2GlobalXY(0),
    h2GlobalSmearedXY(0),
    h2GlobalDeltaXY(0),
    h3GlobalDeltaXYDisk(0),
    h3GlobalDeltaXYR(0) { }
    
int StFstSlowSimMaker::Init() {

	if(mHist){
		fOut = new TFile(mQAFileName.Data(), "RECREATE");
		AddHist(hTrutHitYXDisk = new TH3F("hTrutHitYXDisk", "Global hits before segmentation", 151, -75.5, 75.5, 151, -75.5, 75.5, 10, 0, 10));
		AddHist(hTrutHitRDisk = new TH2F("hTrutHitRDisk", "Global hits before segmentation", 400, 0, 40, 10, 0, 10));
		AddHist(hTrutHitRShower[0] = new TH2F("hTrutHitRShower_4", "Global hits before segmentation", 400, 0, 40, 20, -10, 10));
		AddHist(hTrutHitRShower[1] = new TH2F("hTrutHitRShower_5", "Global hits before segmentation", 400, 0, 40, 20, -10, 10));
		AddHist(hTrutHitRShower[2] = new TH2F("hTrutHitRShower_6", "Global hits before segmentation", 400, 0, 40, 20, -10, 10));
		AddHist(hMCHit[0] = new TH2F("hMCHitDisk_4", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hMCHit[1] = new TH2F("hMCHitDisk_5", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hMCHit[2] = new TH2F("hMCHitDisk_6", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hMCPhiZOut[0] = new TH2F("hMCPhiZOutDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZOut[1] = new TH2F("hMCPhiZOutDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZOut[2] = new TH2F("hMCPhiZOutDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZIn[0] = new TH2F("hMCPhiZInDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZIn[1] = new TH2F("hMCPhiZInDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZIn[2] = new TH2F("hMCPhiZInDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCHit[0] = new TH2F("hRCHitDisk_4", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hRCHit[1] = new TH2F("hRCHitDisk_5", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hRCHit[2] = new TH2F("hRCHitDisk_6", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hRCPhiZOut[0] = new TH2F("hRCPhiZOutDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZOut[1] = new TH2F("hRCPhiZOutDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZOut[2] = new TH2F("hRCPhiZOutDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZIn[0] = new TH2F("hRCPhiZInDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZIn[1] = new TH2F("hRCPhiZInDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZIn[2] = new TH2F("hRCPhiZInDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hTrutHitPhiDisk = new TH2F("hTrutHitPhiDisk", "Global hits before segmentation", 360, 0, 360, 10, 0, 10));
		AddHist(hTrutHitPhiZ = new TH2F("hTrutHitPhiZ", "Global hits before segmentation", 360, 0, 360, 6000, 0, 600));
		AddHist(hRecoHitYXDisk = new TH3F("hRecoHitYXDisk", "Global hits after segmentation", 151, -75.5, 75.5, 151, -75.5, 75.5, 10, 0, 10));
		AddHist(hRecoHitRDisk = new TH2F("hRecoHitRDisk", "Global hits after segmentation", 400, 0, 40, 10, 0, 10));
		AddHist(hRecoHitPhiDisk = new TH2F("hRecoHitPhiDisk", "Global hits after segmentation", 360, 0, 360, 10, 0, 10));
		AddHist(hRecoHitPhiZ = new TH2F("hRecoHitPhiZ", "Global hits after segmentation", 360, 0, 360, 6000, 0, 600));
		AddHist(hGlobalDRDisk = new TH2F("hGlobalDRDisk", "; Reco. r - MC r [cm]; Events;", 1000, -50, 50, 10, 0, 10));
		AddHist(hGlobalZ = new TH1F("hGlobalZ", "; Z [cm]; Events;", 6000, 0, 600));
		AddHist(h3GlobalDeltaXYR = new TH3F("h3GlobalDeltaXYR", ";globalDeltaX; globalDeltaY; R", 300, -0.3, 0.3, 300, -3, 3, 100, 0, 30));
		AddHist(h2GlobalXY = new TH2F("h2GlobalXY", ";globalX; globalY", 1510, -75.5, 75.5, 1510, -75.5, 75.5));
		AddHist(h2GlobalSmearedXY = new TH2F("h2GlobalSmearedXY", ";globalSmearedX; globalSmearedY", 1510, -75.5, 75.5, 1510, -75.5, 75.5));
		AddHist(h2GlobalDeltaXY = new TH2F("h2GlobalDeltaXY", ";globalDeltaX; globalDeltaY", 151, -75.5, 75.5, 151, -75.5, 75.5));
		AddHist(h3GlobalDeltaXYDisk = new TH3F("h3GlobalDeltaXYDisk", ";globalDeltaX; globalDeltaY; Disk", 151, -75.5, 75.5, 151, -75.5, 75.5, 10, 0, 10));
	        for ( int id = 0; id < 3; id++ ){ 
                        std::string HistName = Form("hMCHitCount_Disk%d_InnerSi", id);
                        AddHist(hMCHitCountDisk[id][0] = new TH2F(HistName.c_str(), HistName.c_str(), 120, 4.0, 17.5, 180, -0.1, 2*TMath::Pi()+0.1));
                        HistName = Form("hRCHitCount_Disk%d_InnerSi", id);
                        AddHist(hRCHitCountDisk[id][0] = new TH2F(HistName.c_str(), HistName.c_str(), 120, 4.0, 17.5, 180, -0.1, 2*TMath::Pi()+0.1));
                        HistName = Form("hMCHitCount_Disk%d_OuterSi", id);
                        AddHist(hMCHitCountDisk[id][1] = new TH2F(HistName.c_str(), HistName.c_str(), 120, 15.5, 29.0, 180, -0.1, 2*TMath::Pi()+0.1));
                        HistName = Form("hRCHitCount_Disk%d_OuterSi", id);
                        AddHist(hRCHitCountDisk[id][1] = new TH2F(HistName.c_str(), HistName.c_str(), 120, 15.5, 29.0, 180, -0.1, 2*TMath::Pi()+0.1));
                        for ( int ir = 0; ir < 8; ir++ ){
                                HistName = Form("hRResiduals_Disk%d_Rstrip%d", id, ir);
                                AddHist(hRResiduals[id][ir] = new TH1F(HistName.c_str(), HistName.c_str(), 50, -16.0, 16.0));
                                HistName = Form("hPhiResiduals_Disk%d_Rstrip%d", id, ir);
                                AddHist(hPhiResiduals[id][ir] = new TH1F(HistName.c_str(), HistName.c_str(), 200, -0.01, 0.01)); 
                        }
                }
        }
	return StMaker::Init();
}

void StFstSlowSimMaker::SetDisk(const int i, const float rmn, const float rmx) {
    FstGlobalSlow::RMIN[i] = rmn;
    FstGlobalSlow::RMAX[i] = rmx;
}

Int_t StFstSlowSimMaker::Make() {
	LOG_INFO << "StFstSlowSimMaker::Make" << endm;
        
	// Get the existing StEvent, or add one if it doesn't exist.
	StEvent *event = static_cast<StEvent *>(GetDataSet("StEvent"));
	if (!event) {
		event = new StEvent;
		AddData(event);
		LOG_DEBUG << "Creating StEvent" << endm;
	}

	if (0 == event->rndHitCollection()) {
		event->setRnDHitCollection(new StRnDHitCollection());
		LOG_DEBUG << "Creating StRnDHitCollection for FTS" << endm;
	}

	// Digitize GEANT FTS hits
	FillSilicon(event);

	return kStOk;
}

/* Fill an event with StFtsHits. */
/* This should fill StFtsStrip for realistic simulator and let clustering fill StFtsHit */
/* For now skipping StFtsStrip and clustering, and fill StFtsHits directly here*/

void StFstSlowSimMaker::FillSilicon(StEvent *event) {

	StRnDHitCollection *fsicollection = event->rndHitCollection();

	const int NDISC = 6;
	const int MAXR = mNumR;
	const int MAXPHI = mNumPHI * mNumSEC;

	float X0[] = {0, 0, 0, 0, 0, 0};
	float Y0[] = {0, 0, 0, 0, 0, 0};
	
	if (mRaster > 0)
		for (int i = 0; i < 6; i++) {
			X0[i] = mRaster * TMath::Cos(i * 60 * TMath::DegToRad());
			Y0[i] = mRaster * TMath::Sin(i * 60 * TMath::DegToRad());
		}

	
	// maps for hit and energy for each disk's r-phi strip
	std::map< FstGlobalSlow::FstKeyTriple, StRnDHit* > hitMap;
	std::map< FstGlobalSlow::FstKeyTriple, double > energySum;
	std::map< FstGlobalSlow::FstKeyTriple, double > energyMax;

	// Read the g2t table
	St_g2t_fts_hit *hitTable = static_cast<St_g2t_fts_hit *>(GetDataSet("g2t_fsi_hit"));
	if (!hitTable) {
		LOG_INFO << "g2t_fsi_hit table is empty" << endm;
		return; // Nothing to do
	}

	const Int_t nHits = hitTable->GetNRows();
	LOG_DEBUG << "g2t_fsi_hit table has " << nHits << " hits" << endm;
	const g2t_fts_hit_st *hit = hitTable->GetTable();
	
	StPtrVecRnDHit hits;
        
        std::vector<TVector2> mchits;

	// track table
	St_g2t_track *trkTable = static_cast<St_g2t_track *>(GetDataSet("g2t_track"));
	if (!trkTable) {
		LOG_INFO << "g2t_track table is empty" << endm;
		return; // Nothing to do
	}

	const Int_t nTrks = trkTable->GetNRows();
	
	LOG_DEBUG << "g2t_track table has " << nTrks << " tracks" << endm;
	
	const g2t_track_st *trk = trkTable->GetTable();
	
	int count = 0;
	for (Int_t i = 0; i < nHits; ++i) {
		
		hit = (g2t_fts_hit_st *)hitTable->At(i);
		
		// skip bad hits
		if (!hit) 
			continue;


		int volume_id = hit->volume_id;
		int disk = volume_id / 1000;         // disk id
		int wedge = (volume_id % 1000) / 10; // wedge id
		int sensor = volume_id % 10;         // sensor id
		
		// used as an index for various arrays
		size_t disk_index = disk - 1;
                if (FstGlobalSlow::verbose){
		    LOG_INFO << "disk = " << disk << ", wedge = " << wedge << ", sensor = " << sensor << endm;
                }
		// skip non-FST hits
		if ( disk > 6 ) continue;
                
		double energy = hit->de;
		int track = hit->track_p;

		trk = (g2t_track_st *)trkTable->At(track);
		int isShower = false;
		if (trk)
			isShower = trk->is_shower;

		// raster coordinate offsets
		double xc = X0[disk_index];
		double yc = Y0[disk_index];

		// hit coordinates
		double x = hit->x[0];
		double y = hit->x[1];
		double z = hit->x[2];

		if (z > 200)
			continue; // skip large disks

		// rastered
		double rastered_x = x - xc;
		double rastered_y = y - yc;

		double r = sqrt(x * x + y * y);
		double p = atan2(y, x);

		// rastered
		double rr = sqrt(rastered_x * rastered_x + rastered_y * rastered_y);
		double pp = atan2(rastered_y, rastered_x);


		// wrap an angle between 0 and 2pi
		auto wrapAngle = [&]( double angle ) {
			angle = fmod( angle, 2.0 * M_PI );
			if ( angle < 0 )
				angle += 2.0 * M_PI;
			return angle;
		};

		p = wrapAngle( p );
		pp = wrapAngle( pp );

		LOG_DEBUG << "rr = " << rr << " pp=" << pp << endm;
		LOG_DEBUG << "RMIN = " << FstGlobalSlow::RMIN[disk_index] << " RMAX= " << FstGlobalSlow::RMAX[disk_index] << endm;

		// Cuts made on rastered value to require the r value is within limits
		if (rr < FstGlobalSlow::RMIN[disk_index] || rr > FstGlobalSlow::RMAX[disk_index])
			continue;

		LOG_DEBUG << "rr = " << rr << endm;

		// Strip numbers on rastered value
		int r_index = floor(MAXR * (rr - FstGlobalSlow::RMIN[disk_index]) / (FstGlobalSlow::RMAX[disk_index] - FstGlobalSlow::RMIN[disk_index]));
		
		// this gives a different conflicting answer for r_index and does not handle r outside of range
		for (int ii = 0; ii < MAXR; ii++)
			if (rr > FstGlobalSlow::RSegment[ii] && rr <= FstGlobalSlow::RSegment[ii + 1])
				r_index = ii;
 	
		// Phi number
		int phi_index = int(MAXPHI * pp / 2.0 / M_PI);

		if (r_index >= 8)
			continue;
                
                int mcr_index = r_index;

                // uses measured cross talk values on FST Modules (UIC Cosmic Ray Test Stand) to shift the r_index of the rastered hit
                GetCrossTalkRSegment( r_index ) ;

		if (MAXR)
			assert(r_index < MAXR);
		if (MAXPHI)
			assert(phi_index < MAXPHI);

		StRnDHit *fsihit = nullptr;

		// key of this disk's r & phi strip
		auto threeKey = std::tie( disk_index, r_index, phi_index );
		
		if (hitMap.count( threeKey ) == 0) { // New hit

			if (FstGlobalSlow::verbose){
				LOG_INFO << Form("NEW d=%1d xyz=%8.4f %8.4f %8.4f r=%8.4f phi=%8.4f iR=%2d iPhi=%4d dE=%8.4f[MeV] truth=%d",
						disk, x, y, z, r, p, r_index, phi_index, energy * 1000.0, track)
					<< endm;
			}

			count++;
			fsihit = new StRnDHit();
			fsihit->setDetectorId(kFtsId);
			fsihit->setLayer(disk);
                        fsihit->setLadder(wedge-1);
			fsihit->setWafer(sensor-1);
                        fsihit->setVolumeId(volume_id); // volume id of MC hit from Geant
                        fsihit->setExtraByte0(mcr_index); // MC radial strip index
                        fsihit->setExtraByte1(r_index); // RC radial strip index after cross talk

			//
			// Set position and position error based on radius-constant bins
			//
			double p0 = (phi_index + 0.5) * 2.0 * M_PI / double(MAXPHI);
			double dp = 2.0 * M_PI / double(MAXPHI) / FstGlobalSlow::SQRT12;

			// ONLY valid for the disk array 456, no difference for each disk
			double r0 = (FstGlobalSlow::RSegment[r_index] + FstGlobalSlow::RSegment[r_index + 1]) * 0.5;
			double dr = FstGlobalSlow::RSegment[r_index + 1] - FstGlobalSlow::RSegment[r_index];

                        //LOG_INFO << "NEW R: " << r0 << "     NEW PHI: " << p0 << endm;
			
			double x0 = r0 * cos(p0) + xc;
			double y0 = r0 * sin(p0) + yc;
			assert(TMath::Abs(x0) + TMath::Abs(y0) > 0);
			double dz = 0.03 / FstGlobalSlow::SQRT12;
			double er = dr / FstGlobalSlow::SQRT12;
			fsihit->setPosition(StThreeVectorF(x0, y0, z));
			
			fsihit->setPositionError(StThreeVectorF(er, dp, dz));
			// set covariance matrix
			fsihit->setErrorMatrix(&FstGlobalSlow::Hack1to6(fsihit)[0][0]);

			fsihit->setCharge(energy);
			fsihit->setIdTruth(track, 100);
			hits.push_back(fsihit);
			
			hitMap[ threeKey ] = fsihit;
			energySum[ threeKey ] = energy;
			energyMax[ threeKey ] = energy;

			if (FstGlobalSlow::verbose){
				LOG_INFO << Form("NEW d=%1d xyz=%8.4f %8.4f %8.4f ", disk, x, y, z) << endm;
				LOG_INFO << Form("smeared xyz=%8.4f %8.4f %8.4f ", fsihit->position().x(), fsihit->position().y(), fsihit->position().z()) << endm;
			}

			if(mHist){
				TVector2 hitpos_mc(x, y);
                                mchits.push_back(hitpos_mc);
				TVector2 hitpos_rc(fsihit->position().x(), fsihit->position().y());

				hTrutHitYXDisk->Fill(x, y, disk);
				hTrutHitRDisk->Fill(hitpos_mc.Mod(), disk);
			
                                double rrc = hitpos_rc.Mod();
                                double prc = hitpos_rc.Phi();
                                double rmc = hitpos_mc.Mod();
                                double pmc = hitpos_mc.Phi();
	                        
				if (disk == 4)
					hTrutHitRShower[0]->Fill(hitpos_mc.Mod(), isShower);
				if (disk == 5)
					hTrutHitRShower[1]->Fill(hitpos_mc.Mod(), isShower);
				if (disk == 6)
					hTrutHitRShower[2]->Fill(hitpos_mc.Mod(), isShower);

                                
                                hMCHit[disk-4]->Fill(x,y);
                                hRCHit[disk-4]->Fill(fsihit->position().x(), fsihit->position().y());

                                if (sensor == 1){
                                  hMCPhiZIn[disk-4]->Fill(pmc * 180.0 / TMath::Pi(), z);
                                  hRCPhiZIn[disk-4]->Fill(prc * 180.0 / TMath::Pi(), z); 
                                }
                                else {
                                  hMCPhiZOut[disk-4]->Fill(pmc * 180.0 / TMath::Pi(), z);
                                  hRCPhiZOut[disk-4]->Fill(prc * 180.0 / TMath::Pi(), z);
                                }

				hTrutHitPhiDisk->Fill(hitpos_mc.Phi() * 180.0 / TMath::Pi(), disk);
				hTrutHitPhiZ->Fill(hitpos_mc.Phi() * 180.0 / TMath::Pi(), z);
				hRecoHitYXDisk->Fill(fsihit->position().x(), fsihit->position().y(), disk);
				hRecoHitRDisk->Fill(hitpos_rc.Mod(), disk);
				hRecoHitPhiDisk->Fill(hitpos_rc.Phi() * 180.0 / TMath::Pi(), disk);
				hRecoHitPhiZ->Fill(hitpos_rc.Phi() * 180.0 / TMath::Pi(), z);
				hGlobalDRDisk->Fill(hitpos_rc.Mod() - hitpos_mc.Mod(), disk);
				hGlobalZ->Fill(fsihit->position().z());

				//cout << "CHECK : " << fsihit->position().x()-x << " |||  "<<  fsihit->position().y()-y << endl;
				h2GlobalXY->Fill(x, y);
				h2GlobalSmearedXY->Fill(fsihit->position().x(), fsihit->position().y());
				h2GlobalDeltaXY->Fill(fsihit->position().x() - x, fsihit->position().y() - y);
				h3GlobalDeltaXYDisk->Fill(fsihit->position().x() - x, fsihit->position().y() - y, disk);

				h3GlobalDeltaXYR->Fill(fsihit->position().x() - x, fsihit->position().y() - y, sqrt(pow(fsihit->position().x(), 2) + pow(fsihit->position().y(), 2)));
			        hRResiduals[disk][r_index]->Fill(hitpos_rc.Mod() - hitpos_mc.Mod());
				hPhiResiduals[disk][r_index]->Fill(hitpos_rc.Phi() - hitpos_mc.Phi());    
                        }
		}
		else { // Hit on this strip already exists, adding energy to old hit
			// get hit from the map
			fsihit = hitMap[ threeKey ] ;
			fsihit->setCharge(fsihit->charge() + energy);

			// Add energy to running sum
			energySum[ threeKey ] = energySum[ threeKey ] + energy;

			if (energy> energyMax[ threeKey ])
				energyMax[ threeKey ] = energy;

			// keep idtruth but dilute it...
			track = fsihit->idTruth();

			fsihit->setIdTruth(track, 100 * (energyMax[ threeKey ] / energySum[ threeKey ]));
		}
	} // loop on hits
	int nfsihit = hits.size();

	StarRandom &rand = StarRandom::Instance();
         
        auto wrapAngle = [&]( double angle ) {
                        angle = fmod( angle, 2.0 * M_PI );
                        if ( angle < 0 )
                                angle += 2.0 * M_PI;
                        return angle;
                };

	// skip adding hits based on Efficiency of Inner or Outer silicon sensors (inner: sensor = 1, outer: sensor = 2 || 3)
        for (int i = 0; i < nfsihit; i++) { 
                int sensor = hits[i]->volumeId() % 10; // sensor id
                int disk = (hits[i]->volumeId() / 1000) - 4;
                int sensor_index = 0;
                if (sensor > 1) sensor_index = 1;
 
                if ( mHist ){
                        hMCHitCountDisk[disk][sensor_index]->Fill( mchits[i].Mod(), mchits[i].Phi() );
                }

                bool save_hit = true;
                if ( mInEff ){
                        double rand_save = rand.flat();
                        if ( (sensor == 1 && rand_save < mInnerSiInEff) || (sensor > 1 && rand_save < mOuterSiInEff) )
                                save_hit = false; 
                }
                if ( mHist && save_hit ) {
                        hRCHitCountDisk[disk][sensor_index]->Fill( mchits[i].Mod(), mchits[i].Phi() );
                        hRResiduals[disk][hits[i]->extraByte0()]->Fill( sqrt(pow(hits[i]->position().x(), 2) + pow(hits[i]->position().y(), 2)) - mchits[i].Mod() );
                        hPhiResiduals[disk][hits[i]->extraByte0()]->Fill( wrapAngle(hits[i]->position().phi()) -  mchits[i].Phi() );
                        //LOG_INFO << "phi residual: " << wrapAngle(hits[i]->position().phi()) -  mchits[i].Phi() << endm;
                }
                if ( save_hit ) fsicollection->addHit(hits[i]);	
	}
	if (FstGlobalSlow::verbose) {
		LOG_DEBUG << Form("Found %d/%d g2t hits in %d cells, created %d hits with ADC>0", count, nHits, nfsihit, fsicollection->numberOfHits()) << endm;
	}

}
//

int StFstSlowSimMaker::Finish() {
	if(mHist){
		fOut->cd();
		hTrutHitYXDisk->Write();
		hTrutHitRDisk->Write();
		hTrutHitRShower[0]->Write();
		hTrutHitRShower[1]->Write();
		hTrutHitRShower[2]->Write();
		hTrutHitPhiDisk->Write();
		hTrutHitPhiZ->Write();
		hRecoHitYXDisk->Write();
		hRecoHitRDisk->Write();
		hRecoHitPhiDisk->Write();
		hRecoHitPhiZ->Write();
		hGlobalDRDisk->Write();
		hGlobalZ->Write();
		h3GlobalDeltaXYR->Write();
		h2GlobalXY->Write();
		h2GlobalSmearedXY->Write();
		h2GlobalDeltaXY->Write();
		h3GlobalDeltaXYDisk->Write();
                hMCHit[0]->Write();
                hMCHit[1]->Write();
                hMCHit[2]->Write();
                hMCPhiZOut[0]->Write();
                hMCPhiZOut[1]->Write();
                hMCPhiZOut[2]->Write();
                hMCPhiZIn[0]->Write();
                hMCPhiZIn[1]->Write();
                hMCPhiZIn[2]->Write();
                hRCHit[0]->Write();
                hRCHit[1]->Write();
                hRCHit[2]->Write();
                hRCPhiZOut[0]->Write();
                hRCPhiZOut[1]->Write();
                hRCPhiZOut[2]->Write();
                hRCPhiZIn[0]->Write();
                hRCPhiZIn[1]->Write();
                hRCPhiZIn[2]->Write();
              
                for ( int id = 0; id < 3; id++ ){
                        for ( int is = 0; is < 2; is++ ){
                                hMCHitCountDisk[id][is]->Write();
                                hRCHitCountDisk[id][is]->Write();
                        }
                        for ( int ir = 0; ir < 8; ir++ ){
                                hRResiduals[id][ir]->Write();
                                hPhiResiduals[id][ir]->Write();
                        }
                } 
		fOut->Close();
	}
	return kStOK;
}

void StFstSlowSimMaker::GetCrossTalkRSegment( int &rstrip ){

  	const double crossTalkRate[8][7] = {
        	{0.0000, 0.0000, 0.0000, 0.7419, 0.9413, 0.9775, 1.0000},
        	{0.0000, 0.0000, 0.0247, 0.8773, 0.9809, 1.0000, 1.0000},
        	{0.0000, 0.0052, 0.0668, 0.8726, 1.0000, 1.0000, 1.0000},
        	{0.0097, 0.0651, 0.1495, 1.0000, 1.0000, 1.0000, 1.0000},
        	{0.0000, 0.0000, 0.0000, 0.8780, 0.9810, 0.9915, 1.0000},
        	{0.0000, 0.0000, 0.1343, 0.9141, 0.9889, 1.0000, 1.0000},
        	{0.0000, 0.0372, 0.1889, 0.9288, 1.0000, 1.0000, 1.0000},
                {0.0415, 0.2166, 0.3871, 1.0000, 1.0000, 1.0000, 1.0000}
        };

        StarRandom &rand = StarRandom::Instance();
        double rand_ct = rand.flat();        

  	for (int i_delta = 0; i_delta < 7; ++i_delta){
      		if (rand_ct >= crossTalkRate[rstrip][i_delta] && rand_ct < crossTalkRate[rstrip][i_delta+1]){
        		rstrip += ( i_delta - 2 );
                        return;
      		}
    	}
    	if (rand_ct >= 0.0 && rand_ct < crossTalkRate[rstrip][0]){
      		rstrip -= 3;
                return;
    	}
}
  
//_____________________________________________________________________________
StMatrixF FstGlobalSlow::Hack1to6(const StHit *stHit) {
	//   X = R*cos(Fi), Y=R*sin(Fi), Z = z
	//   dX/dR  = (    cos(Fi)  ,sin(Fi),0)
	//   dX/dFi = (-R*sin(Fi), R*cos(Fi),0)
	//   dX/dZ  = (         0,         0,1)

	auto hiPos = stHit->position();
	auto hiErr = stHit->positionError();
	double Rxy = sqrt(hiPos[0] * hiPos[0] + hiPos[1] * hiPos[1]);
	double cosFi = hiPos[0] / Rxy;
	double sinFi = hiPos[1] / Rxy;
	double T[3][3] = {{cosFi, -Rxy * sinFi, 0}, {sinFi, Rxy * cosFi, 0}, {0, 0, 1}};
	double Ginp[6] = {hiErr[0] * hiErr[0], 0, hiErr[1] * hiErr[1], 0, 0, hiErr[2] * hiErr[2]};
	double Gout[6];

	TCL::trasat(T[0], Ginp, Gout, 3, 3);
	StMatrixF mtxF(3, 3);

	for (int i = 0, li = 0; i < 3; li += ++i) {
		for (int j = 0; j <= i; j++) {
			mtxF[i][j] = Gout[li + j];
			mtxF[j][i] = mtxF[i][j];
		}
	}

	return mtxF;
}
