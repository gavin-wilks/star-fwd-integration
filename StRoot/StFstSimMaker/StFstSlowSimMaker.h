#ifndef ST_FST_SLOW_SIM_MAKER_H
#define ST_FST_SLOW_SIM_MAKER_H

class StFtsHit;
class StEvent;
class StRnDHitCollection;
class StRnDHit;

#include "StChain/StMaker.h"
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


class StFstSlowSimMaker : public StMaker {
	public:
		explicit StFstSlowSimMaker(const Char_t *name = "fstFastSim");
		virtual ~StFstSlowSimMaker() {}
		int Make();
		int Init();
		int Finish();
		virtual const char *GetCVS() const;
		
		/// Set offset for each disk ( x=R*cos(idisk*60 degrees), y=R*sin(...) )
		void SetRaster(float R = 1.0) { mRaster = R; }

		/// Set min/max active radii for each disk
		void SetDisk(const int i, const float rmn, const float rmx);
		void SetInEfficiency(bool ineff = false, float innersi = 0.0, float outersi = 0.0) { mInEff = ineff; mInnerSiInEff = innersi; mOuterSiInEff = outersi; }
		void SetQAFileName(TString filename = "") { mQAFileName = filename; }
		void SetFillHist(bool hist = false) { mHist = hist; }
                void SetCrossTalk( bool crosstalk = false) { mCrossTalk = crosstalk; }		

	private:
		void FillSilicon(StEvent *event);
		StRnDHitCollection *hitCollection = nullptr;

                void GetCrossTalkRSegment(int &rstrip);             

		int mNumR;
		int mNumPHI;
		int mNumSEC;
		float mRaster;
		bool mInEff;
                float mInnerSiInEff;
                float mOuterSiInEff;
                bool mCrossTalk;
		bool mHist;
		TString mQAFileName;

		TH3F *hTrutHitYXDisk;
		TH2F *hTrutHitRDisk;
		TH2F *hTrutHitRShower[3];
                TH2F *hMCHit[3];
                TH2F *hMCPhiZOut[3];
                TH2F *hMCPhiZIn[3];
                TH2F *hRCHit[3];
                TH2F *hRCPhiZOut[3];
                TH2F *hRCPhiZIn[3];

		TH2F *hTrutHitPhiDisk;
		TH2F *hTrutHitPhiZ;
		TH3F *hRecoHitYXDisk;
		TH2F *hRecoHitRDisk;
		TH2F *hRecoHitPhiDisk;
		TH2F *hRecoHitPhiZ;
		TH2F *hGlobalDRDisk;
		TH1F *hGlobalZ;

		TH2F *h2GlobalXY;
		TH2F *h2GlobalSmearedXY;
		TH2F *h2GlobalDeltaXY;
		TH3F *h3GlobalDeltaXYDisk;
		TH3F *h3GlobalDeltaXYR;
                TH2F *hMCHitCountDisk[3][2];
                TH2F *hRCHitCountDisk[3][2];
                TH1F *hRResiduals[3][8];
                TH1F *hPhiResiduals[3][8];
		TFile *fOut;

		ClassDef(StFstSlowSimMaker, 0)
};

inline const char *StFstSlowSimMaker::GetCVS() const {
	static const char cvs[] = "Tag $Name:  $ $Id: StFstSlowSimMaker.h,v 1.1 2018/11/06 18:56:05 jdb Exp $ built " __DATE__ " " __TIME__;
	return cvs;
}

#endif
