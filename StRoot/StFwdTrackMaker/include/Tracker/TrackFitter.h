#ifndef TRACK_FITTER_H
#define TRACK_FITTER_H

#include "GenFit/ConstField.h"
#include "GenFit/EventDisplay.h"
#include "GenFit/Exception.h"
#include "GenFit/FieldManager.h"
#include "GenFit/KalmanFitStatus.h"
#include "GenFit/KalmanFitter.h"
#include "GenFit/KalmanFitterInfo.h"
#include "GenFit/KalmanFitterRefTrack.h"
#include "GenFit/MaterialEffects.h"
#include "GenFit/PlanarMeasurement.h"
#include "GenFit/RKTrackRep.h"
#include "GenFit/SpacepointMeasurement.h"
#include "GenFit/StateOnPlane.h"
#include "GenFit/TGeoMaterialInterface.h"
#include "GenFit/Track.h"
#include "GenFit/TrackPoint.h"

#include "Criteria/SimpleCircle.h"

#include "TDatabasePDG.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"

#include <vector>
#include <memory>

#include "StFwdTrackMaker/Common.h"

#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/STARField.h"
#include "StFwdTrackMaker/include/Tracker/loguru.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"

#include "StarGenerator/UTIL/StarRandom.h"

/* Cass for fitting tracks(space points) with GenFit
 *
 */
class TrackFitter {

  public:
    // ctor 
    // provide the main configuration object
    TrackFitter(FwdTrackerConfig _mConfig) : mConfig(_mConfig) {
        mTrackRep = 0;
        mFitTrack = 0;
    }


    void setup() {

        // the geometry manager that GenFit will use
        TGeoManager * gMan = nullptr;

        // Setup the Geometry used by GENFIT
        TGeoManager::Import(mConfig.get<string>("Geometry", "fGeom.root").c_str());
        gMan = gGeoManager;
        // Set up the material interface and set material effects on/off from the config
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        genfit::MaterialEffects::getInstance()->setNoEffects(mConfig.get<bool>("TrackFitter::noMaterialEffects", false)); // false means defaul ON

        // Determine which Magnetic field to use
        // Either constant field or real field from StarFieldAdaptor
        if (mConfig.get<bool>("TrackFitter:constB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 5.)); // 0.5 T Bz
            LOG_INFO << "StFwdTrackMaker: Tracking with constant magnetic field" << endl;
        } else {
            mBField = std::unique_ptr<genfit::AbsBField>(new StarFieldAdaptor());
            LOG_INFO << "StFwdTrackMaker: Tracking with StarFieldAdapter" << endl;
        }
        // we must have one of the two available fields at this point
        genfit::FieldManager::getInstance()->init(mBField); 

        // initialize the main mFitter using a KalmanFitter with reference tracks
        mFitter = std::unique_ptr<genfit::AbsKalmanFitter>(new genfit::KalmanFitterRefTrack());

        // Here we load several options from the config, 
        // to customize the mFitter behavior
        mFitter->setMaxFailedHits(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxFailedHits", -1)); // default -1, no limit
        mFitter->setDebugLvl(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:DebugLvl", 0)); // default 0, no output
        mFitter->setMaxIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxIterations", 4)); // default 4 iterations
        mFitter->setMinIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MinIterations", 0)); // default 0 iterations

        // FwdGeomUtils looks into the loaded geometry and gets detector z locations if present
        FwdGeomUtils fwdGeoUtils( gMan );

        // these default values are the default if the detector is 
        // a) not found in the geometry 
        // b) not provided in config

        // NOTE: these defaults are needed since the geometry file might not include FST (bug being worked on separately)
        mFSTZLocations = mConfig.getVector<float>("TrackFitter.Geometry:fst", {140.286011, 154.286011, 168.286011 });
        if ( fwdGeoUtils.siZ( 0 ) > 1.0 ) { // returns 0.0 on failure
            // The geometry has the silicon detector, so use its z locations
            mFSTZLocations.clear();
            mFSTZLocations.push_back( fwdGeoUtils.siZ( 0 ) );
            mFSTZLocations.push_back( fwdGeoUtils.siZ( 1 ) );
            mFSTZLocations.push_back( fwdGeoUtils.siZ( 2 ) );
        } else {
            LOG_WARNING << "Using FST z-locations from config or defautl, may not match hits" << endm;
        }

        // Now add the Si detector planes at the desired location
        for (auto z : mFSTZLocations) {
            mFSTPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0, 0, z), TVector3(1, 0, 0), TVector3(0, 1, 0) )
                )
            );
        }

        // Now load FTT
        // mConfig.getVector<>(...) requires a default, hence the 
        mFTTZLocations = mConfig.getVector<float>("TrackFitter.Geometry:ftt", {0.0f, 0.0f, 0.0f, 0.0f});

        if ( fwdGeoUtils.stgcZ( 0 ) > 1.0 ) { // returns 0.0 on failure
            mFTTZLocations.clear();

            mFTTZLocations.push_back( fwdGeoUtils.stgcZ( 0 ) );
            mFTTZLocations.push_back( fwdGeoUtils.stgcZ( 1 ) );
            mFTTZLocations.push_back( fwdGeoUtils.stgcZ( 2 ) );
            mFTTZLocations.push_back( fwdGeoUtils.stgcZ( 3 ) );
        } else {
            LOG_WARNING << "Using FTT z-locations from config or default, may not match hits" << endm;
        }

        if ( mFTTZLocations.size() != 4 ){
            LOG_ERROR << "Wrong number of FTT layers, got " << mFTTZLocations.size() << " but expected 4" << endm;
        }

        for (auto z : mFTTZLocations) {
            mFTTPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0, 0, z), TVector3(1, 0, 0), TVector3(0, 1, 0))
                )
            );
        }

        // get default vertex values used in simulation from the config
        mVertexSigmaXY = mConfig.get<float>("TrackFitter.Vertex:sigmaXY", 1.0f);
        mVertexSigmaZ = mConfig.get<float>("TrackFitter.Vertex:sigmaZ", 30.0f);
        mVertexPos = mConfig.getVector<float>("TrackFitter.Vertex:pos", {0.0f,0.0f,0.0f});
        mIncludeVertexInFit = mConfig.get<bool>("TrackFitter.Vertex:includeInFit", false);

        if ( mGenHistograms )
            makeHistograms();
    }

    void makeHistograms() {
        std::string n = "";
        mHist["ECalProjPosXY"] = new TH2F("ECalProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["ECalProjSigmaXY"] = new TH2F("ECalProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 50, 0, 0.5, 50, 0, 0.5);
        mHist["ECalProjSigmaR"] = new TH1F("ECalProjSigmaR", ";#sigma_{XY} (cm) at ECAL", 50, 0, 0.5);

        mHist["SiProjPosXY"] = new TH2F("SiProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["SiProjSigmaXY"] = new TH2F("SiProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 150, 0, 15, 150, 0, 15);

        mHist["VertexProjPosXY"] = new TH2F("VertexProjPosXY", ";X;Y", 100, -5, 5, 100, -5, 5);
        mHist["VertexProjSigmaXY"] = new TH2F("VertexProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 150, 0, 20, 150, 0, 20);

        mHist["VertexProjPosZ"] = new TH1F("VertexProjPosZ", ";Z;", 100, -50, 50);
        mHist["VertexProjSigmaZ"] = new TH1F("VertexProjSigmaZ", ";#sigma_{Z};", 100, 0, 10);

        mHist["SiWrongProjPosXY"] = new TH2F("SiWrongProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["SiWrongProjSigmaXY"] = new TH2F("SiWrongProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 50, 0, 0.5, 50, 0, 0.5);

        mHist["SiDeltaProjPosXY"] = new TH2F("SiDeltaProjPosXY", ";X;Y", 1000, 0, 20, 1000, 0, 20);

        n = "seed_curv";
        mHist[n] = new TH1F(n.c_str(), ";curv", 1000, 0, 10000);
        n = "seed_pT";
        mHist[n] = new TH1F(n.c_str(), ";pT (GeV/c)", 500, 0, 10);
        n = "seed_eta";
        mHist[n] = new TH1F(n.c_str(), ";eta", 500, 0, 5);

        n = "delta_fit_seed_pT";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) pT (GeV/c)", 500, -5, 5);
        n = "delta_fit_seed_eta";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) eta", 500, 0, 5);
        n = "delta_fit_seed_phi";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) phi", 500, -5, 5);

        n = "FitStatus";
        mHist[n] = new TH1F(n.c_str(), ";", 5, 0, 5);
        FwdTrackerUtils::labelAxis(mHist[n]->GetXaxis(), {"Total", "Pass", "Fail", "GoodCardinal", "Exception"});

        n = "FitDuration";
        mHist[n] = new TH1F(n.c_str(), "; Duraton (ms)", 5000, 0, 50000);

        n = "FailedFitDuration";
        mHist[n] = new TH1F(n.c_str(), "; Duraton (ms)", 500, 0, 50000);
    }

    // writes mHistograms stored in map only if mGenHistograms is true
    void writeHistograms() {
        if ( !mGenHistograms )
            return;
        for (auto nh : mHist) {
            nh.second->SetDirectory(gDirectory);
            nh.second->Write();
        }
    }

    /* Convert the 3x3 covmat to 2x2 by dropping z
    *
    */
    TMatrixDSym CovMatPlane(KiTrack::IHit *h){
        TMatrixDSym cm(2);
        cm(0, 0) = static_cast<FwdHit*>(h)->_covmat(0, 0);
        cm(1, 1) = static_cast<FwdHit*>(h)->_covmat(1, 1);
        cm(0, 1) = static_cast<FwdHit*>(h)->_covmat(0, 1);
        return cm;
    }


    /* FitSimpleCircle
     * Used to determine a seed transverse momentum based on space points
     * Takes a list of space points KiTrack::IHit *
     * Takes three indecise used to lookup three of the possible hits within the list
     */ 
    float fitSimpleCircle(vector<KiTrack::IHit *> trackCand, size_t i0, size_t i1, size_t i2) {
        float curv = 0;

        // ensure that no index is outside of range for FST or FTT volumes
        if (i0 > 12 || i1 > 12 || i2 > 12)
            return 0;

        try {
            KiTrack::SimpleCircle sc(trackCand[i0]->getX(), trackCand[i0]->getY(), trackCand[i1]->getX(), trackCand[i1]->getY(), trackCand[i2]->getX(), trackCand[i2]->getY());
            curv = sc.getRadius();
        } catch (KiTrack::InvalidParameter &e) {
            // if we got here we failed to get  a valid seed. We will still try to move forward but the fit will probably fail
        }

        //  make sure the curv is valid
        if (isinf(curv))
            curv = 999999.9;

        return curv;
    }

    /* seedState
     * Determines the seed position and momentum for a list of space points
     */
    float seedState(vector<KiTrack::IHit *> trackCand, TVector3 &seedPos, TVector3 &seedMom) {
        // we require at least 4 hits,  so this should be gauranteed
        if(trackCand.size() < 3){
            // failure
            return 0.0;
        }
            

        // we want to use the LAST 3 hits, since silicon doesnt have R information
        TVector3 p0, p1, p2;
        // use the closest hit to the interaction point for the seed pos
        FwdHit *hit_closest_to_IP = static_cast<FwdHit *>(trackCand[0]);

        // maps from <key=vol_id> to <value=index in trackCand>
        std::map<size_t, size_t> vol_map; 

        // init the map
        for (size_t i = 0; i < 13; i++)
            vol_map[i] = -1;

        for (size_t i = 0; i < trackCand.size(); i++) {
            auto fwdHit = static_cast<FwdHit *>(trackCand[i]);
            vol_map[abs(fwdHit->_vid)] = i;
            // find the hit closest to IP for the initial position seed
            if (hit_closest_to_IP->getZ() > fwdHit->getZ())
                hit_closest_to_IP = fwdHit;
        }

        // now get an estimate of the pT from several overlapping simple circle fits
        // enumerate the available partitions
        // 12 11 10
        // 12 11 9
        // 12 10 9
        // 11 10 9
        vector<float> curvs;
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[12], vol_map[11], vol_map[10]));
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[12], vol_map[11], vol_map[9]));
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[12], vol_map[10], vol_map[9]));
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[11], vol_map[10], vol_map[9]));

        // average them and exclude failed fits
        float mcurv = 0;
        float nmeas = 0;

        for (size_t i = 0; i < curvs.size(); i++) {
            if (mGenHistograms)
                this->mHist["seed_curv"]->Fill(curvs[i]);
            if (curvs[i] > 10) {
                mcurv += curvs[i];
                nmeas += 1.0;
            }
        }

        if (nmeas >= 1)
            mcurv = mcurv / nmeas;
        else
            mcurv = 10;

        // Now lets get eta information
        // simpler, use farthest points from IP
        if (vol_map[9] < 13)
            p0.SetXYZ(trackCand[vol_map[9]]->getX(), trackCand[vol_map[9]]->getY(), trackCand[vol_map[9]]->getZ());

        if (vol_map[10] < 13)
            p1.SetXYZ(trackCand[vol_map[10]]->getX(), trackCand[vol_map[10]]->getY(), trackCand[vol_map[10]]->getZ());

        const double K = 0.00029979; //K depends on the units used for Bfield
        double pt = mcurv * K * 5; // pT from average measured curv
        double dx = (p1.X() - p0.X());
        double dy = (p1.Y() - p0.Y());
        double dz = (p1.Z() - p0.Z());
        double phi = TMath::ATan2(dy, dx);
        double Rxy = sqrt(dx * dx + dy * dy);
        double theta = TMath::ATan2(Rxy, dz);
        // double eta = -log( tantheta / 2.0 );
        // these starting conditions can probably be improvd, good study for student

        seedMom.SetPtThetaPhi(pt, theta, phi);
        seedPos.SetXYZ(hit_closest_to_IP->getX(), hit_closest_to_IP->getY(), hit_closest_to_IP->getZ());

        if (mGenHistograms) {
            this->mHist["seed_pT"]->Fill(seedMom.Pt());
            this->mHist["seed_eta"]->Fill(seedMom.Eta());
        }

        return mcurv;
    }


    /* RefitTracksWithSiHits
     * Takes a previously fit track re-fits it with the newly added silicon hits 
     * 
     */
    TVector3 refitTrackWithSiHits(genfit::Track *originalTrack, std::vector<KiTrack::IHit *> si_hits) {

        TVector3 pOrig = originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));
        auto cardinalStatus = originalTrack->getFitStatus(originalTrack->getCardinalRep());

        if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
            // in this case the original track did not converge so we should not refit. 
            // probably never get here due to previous checks
            return pOrig;
        }

        // Setup the Track Reps
        auto trackRepPos = new genfit::RKTrackRep(mPdgPiPlus);
        auto trackRepNeg = new genfit::RKTrackRep(mPdgPiMinus);

        // get the space points on the original track
        auto trackPoints = originalTrack->getPointsWithMeasurement();
        
        if ((trackPoints.size() < (mFTTZLocations.size() +1) && mIncludeVertexInFit) || trackPoints.size() < mFTTZLocations.size() ) {
            // we didnt get enough points for a refit
            return pOrig;
        }

        TVectorD rawCoords = trackPoints[0]->getRawMeasurement()->getRawHitCoords();
        double z = mFSTZLocations[0]; //first FTT plane, used if we dont have PV in fit
        if (mIncludeVertexInFit)
            z = rawCoords(2);

        TVector3 seedPos(rawCoords(0), rawCoords(1), z);
        TVector3 seedMom = pOrig;

        // Create the ref track using the seed state
        auto refTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        refTrack->addTrackRep(trackRepNeg);

        genfit::Track &fitTrack = *refTrack;

        size_t firstFTTIndex = 0;
        if (mIncludeVertexInFit) {
            // clone the PRIMARY VERTEX into this track
            fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[0]->getRawMeasurement(), &fitTrack));
            firstFTTIndex = 1; // start on hit index 1 below
        }

        // initialize the hit coords on plane
        TVectorD hitCoords(2);
        hitCoords[0] = 0;
        hitCoords[1] = 0;

        size_t planeId(0);
        int hitId(5);

        // add the hits to the track
        for (auto h : si_hits) {
            if ( nullptr == h ) continue; // if no Si hit in this plane, skip

            hitCoords[0] = h->getX();
            hitCoords[1] = h->getY();
            genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), h->getSector(), ++hitId, nullptr);

            planeId = h->getSector();

            if (mFSTPlanes.size() <= planeId) {
                LOG_WARNING << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
                return pOrig;
            }

            auto plane = mFSTPlanes[planeId];
            measurement->setPlane(plane, planeId);
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
        }
        // start at 0 if PV not included, 1 otherwise 
        for (size_t i = firstFTTIndex; i < trackPoints.size(); i++) {
            // clone the track points into this track
            fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[i]->getRawMeasurement(), &fitTrack));
        }

        try {
            //Track RE-Fit with GENFIT2
            // check consistency of all points
            fitTrack.checkConsistency();

            // do the actual track fit
            mFitter->processTrack(&fitTrack);

            fitTrack.checkConsistency();

            // this chooses the lowest chi2 fit result as cardinal
            fitTrack.determineCardinalRep(); 

        } catch (genfit::Exception &e) {
            // will be caught below by converge check
        }

        if (fitTrack.getFitStatus(fitTrack.getCardinalRep())->isFitConverged() == false) {
            // Did not converge
            return pOrig;
        } else { // we did converge, return new momentum
            
            TVector3 p = fitTrack.getCardinalRep()->getMom(fitTrack.getFittedState(1, fitTrack.getCardinalRep()));
            // get status if needed later
            // auto newStatus = fitTrack.getFitStatus(fitTrack.getCardinalRep());
            return p;
        }
        return pOrig;
    }

    /* Generic method for fitting space points with GenFit
     *
     * 
     */
    TVector3 fitSpacePoints( vector<genfit::SpacepointMeasurement*> spoints, TVector3 &seedPos, TVector3 &seedMom ){
        
        // setup track reps
        auto trackRepPos = new genfit::RKTrackRep(pdg_mu_plus);
        auto trackRepNeg = new genfit::RKTrackRep(pdg_mu_minus);

        // setup track for fit with positive and negative reps
        auto mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack->addTrackRep(trackRepNeg);

        genfit::Track &fitTrack = *mFitTrack;

        // try adding the points to track and fitting
        try {
            for ( size_t i = 0; i < spoints.size(); i++ ){
                fitTrack.insertPoint(new genfit::TrackPoint(spoints[i], &fitTrack));
            }
            // do the fit against the two possible fits
            mFitter->processTrackWithRep(&fitTrack, trackRepPos);
            mFitter->processTrackWithRep(&fitTrack, trackRepNeg);

        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }

        try {
            fitTrack.checkConsistency();

            fitTrack.determineCardinalRep();
            auto cardinalRep = fitTrack.getCardinalRep();

            TVector3 p = cardinalRep->getMom(fitTrack.getFittedState(1, cardinalRep));
            // sucess, return momentum
            return p;
        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }
        return TVector3(0, 0, 0);
    }

    /* Fit a track 
     *
     * 
     */
    TVector3 fitTrack(vector<KiTrack::IHit *> trackCand, double *Vertex = 0, TVector3 *McSeedMom = 0) {
        long long itStart = loguru::now_ns();
        if (mGenHistograms) this->mHist["FitStatus"]->Fill("Total", 1);

        // The PV information, if we want to use it
        TVectorD pv(3);

        if (0 == Vertex) { // randomized from simulation
            pv[0] = mVertexPos[0] + StarRandom::Instance()->gauss(mVertexSigmaXY);
            pv[1] = mVertexPos[1] + StarRandom::Instance()->gauss(mVertexSigmaXY);
            pv[2] = mVertexPos[2] + StarRandom::Instance()->gauss(mVertexSigmaZ);
        } else {
            pv[0] = Vertex[0];
            pv[1] = Vertex[1];
            pv[2] = Vertex[2];
        }

        // get the seed info from our hits
        TVector3 seedMom, seedPos;
        float curv = seedState(trackCand, seedPos, seedMom);

        // create the track representations
        auto trackRepPos = new genfit::RKTrackRep(mPdgPiPlus);
        auto trackRepNeg = new genfit::RKTrackRep(mPdgPiMinus);

        // If we use the PV, use that as the start pos for the track
        if (mIncludeVertexInFit) {
            seedPos.SetXYZ(pv[0], pv[1], pv[2]);
        }

        if (mFitTrack){
            delete mFitTrack;
        }

        if (McSeedMom != nullptr) {
            seedMom = *McSeedMom;
        }

        mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack->addTrackRep(trackRepNeg);

        LOG_INFO
            << "seedPos : (" << seedPos.X() << ", " << seedPos.Y() << ", " << seedPos.Z() << " )"
            << ", seedMom : (" << seedMom.X() << ", " << seedMom.Y() << ", " << seedMom.Z() << " )"
            << ", seedMom : (" << seedMom.Pt() << ", " << seedMom.Eta() << ", " << seedMom.Phi() << " )"
            << endm;

        genfit::Track &fitTrack = *mFitTrack;

        size_t planeId(0);     // detector plane ID
        int hitId(0);       // hit ID

        // initialize the hit coords on plane
        TVectorD hitCoords(2);
        hitCoords[0] = 0;
        hitCoords[1] = 0;

        /******************************************************************************************************************
        * Include the Primary vertex if desired
        ******************************************************************************************************************/
        if (mIncludeVertexInFit) {

            TMatrixDSym hitCov3(3);
            hitCov3(0, 0) = mVertexSigmaXY * mVertexSigmaXY;
            hitCov3(1, 1) = mVertexSigmaXY * mVertexSigmaXY;
            hitCov3(2, 2) = mVertexSigmaZ * mVertexSigmaZ;

            genfit::SpacepointMeasurement *measurement = new genfit::SpacepointMeasurement(pv, hitCov3, 0, ++hitId, nullptr);
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
        }

        /******************************************************************************************************************
		 * loop over the hits, add them to the track
		 ******************************************************************************************************************/
        for (auto h : trackCand) {

            hitCoords[0] = h->getX();
            hitCoords[1] = h->getY();
            
            genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), h->getSector(), ++hitId, nullptr);

            planeId = h->getSector();

            if (mFTTPlanes.size() <= planeId) {
                LOG_WARNING << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
                return TVector3(0, 0, 0);
            }

            auto plane = mFTTPlanes[planeId];
            measurement->setPlane(plane, planeId);
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

            if (abs(h->getZ() - plane->getO().Z()) > 0.05) {
                LOG_WARNING << "Z Mismatch h->z = " << h->getZ() << ", plane->z = "<< plane->getO().Z() <<", diff = " << abs(h->getZ() - plane->getO().Z()) << endm;
            }
        }

        try {
            // do the fit
            mFitter->processTrackWithRep(&fitTrack, trackRepPos);
            mFitter->processTrackWithRep(&fitTrack, trackRepNeg);

        } catch (genfit::Exception &e) {
            if (mGenHistograms) mHist["FitStatus"]->Fill("Exception", 1);
        }

        TVector3 p(0, 0, 0);

        try {
            //check
            fitTrack.checkConsistency();

            // find track rep with smallest chi2
            fitTrack.determineCardinalRep();
            auto cardinalRep = fitTrack.getCardinalRep();
            auto cardinalStatus = fitTrack.getFitStatus(cardinalRep);
            mFitStatus = *cardinalStatus; // save the status of last fit

            // Delete any previous track rep
            if (mTrackRep)
                delete mTrackRep;

            // Clone the cardinal rep for persistency
            mTrackRep = cardinalRep->clone(); // save the result of the fit
            if (fitTrack.getFitStatus(cardinalRep)->isFitConverged() && mGenHistograms ) {
                this->mHist["FitStatus"]->Fill("GoodCardinal", 1);
            }

            if (fitTrack.getFitStatus(trackRepPos)->isFitConverged() == false &&
                fitTrack.getFitStatus(trackRepNeg)->isFitConverged() == false) {
        
                p.SetXYZ(0, 0, 0);
                long long duration = (loguru::now_ns() - itStart) * 1e-6; // milliseconds
                if (mGenHistograms) {
                    this->mHist["FitStatus"]->Fill("Fail", 1);
                    this->mHist["FailedFitDuration"]->Fill(duration);
                }
                return p;
            } // neither track rep converged

            p = cardinalRep->getMom(fitTrack.getFittedState(1, cardinalRep));
            mQ = cardinalRep->getCharge(fitTrack.getFittedState(1, cardinalRep));
            mP = p;

        } catch (genfit::Exception &e) {
            LOG_WARNING << "Exception on track fit: " << e.what() << endm;
            p.SetXYZ(0, 0, 0);

            long long duration = (loguru::now_ns() - itStart) * 1e-6; // milliseconds
            if (mGenHistograms) {
                this->mHist["FitStatus"]->Fill("Exception", 1);
                this->mHist["FailedFitDuration"]->Fill(duration);
            }

            return p;
        } // try/catch 

        long long duration = (loguru::now_ns() - itStart) * 1e-6; // milliseconds
        if (mGenHistograms) {
            this->mHist["FitStatus"]->Fill("Pass", 1);
            this->mHist["delta_fit_seed_pT"]->Fill(p.Pt() - seedMom.Pt());
            this->mHist["delta_fit_seed_eta"]->Fill(p.Eta() - seedMom.Eta());
            this->mHist["delta_fit_seed_phi"]->Fill(p.Phi() - seedMom.Phi());
            this->mHist["FitDuration"]->Fill(duration);
        }
        return p;
    }

    int getCharge() {
        return (int)mQ;
    }

    genfit::FitStatus getStatus() { return mFitStatus; }
    genfit::AbsTrackRep *getTrackRep() { return mTrackRep; }
    genfit::Track *getTrack() { return mFitTrack; }

  protected:
    std::unique_ptr<genfit::AbsBField> mBField;

    FwdTrackerConfig mConfig; // main config object

    // optional histograms, off by default
    std::map<std::string, TH1 *> mHist;
    bool mGenHistograms = false;

    // Main GenFit fitter instance
    std::unique_ptr<genfit::AbsKalmanFitter> mFitter = nullptr;

    // PDG codes for the default plc type for fits
    const int mPdgPiPlus = 211;
    const int mPdgPiMinus = -211;

    // Store the planes for FTT and FST
    vector<genfit::SharedPlanePtr> mFTTPlanes;
    vector<genfit::SharedPlanePtr> mFSTPlanes;

    // det z locations loaded from geom or config
    vector<float> mFSTZLocations, mFTTZLocations;

    // parameter ALIASED from mConfig wrt PV vertex
    float mVertexSigmaXY = 1;
    float mVertexSigmaZ = 30;
    vector<float> mVertexPos;
    bool mIncludeVertexInFit = false;

    // GenFit state
    genfit::FitStatus mFitStatus;
    genfit::AbsTrackRep *mTrackRep;
    genfit::Track *mFitTrack;

    // Fit results
    TVector3 mP;
    double mQ;
};

#endif
