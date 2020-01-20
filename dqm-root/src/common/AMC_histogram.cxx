#include "GEB_histogram.cxx"
#include "TH1.h"

//!A class creating histograms for AMC data
class AMC_histogram: public Hardware_histogram
{
public:
  //!Constructor calls the base constructor of Hardware_histogram. Requires a string for filename, directory, and another string.
  AMC_histogram(const std::string & filename, TDirectory * dir, const std::string & hwid):Hardware_histogram(filename, dir, hwid){}//call base constructor
  ~AMC_histogram()
  {
    delete[] m_gebsH;
  }
  
  //!Books histograms for AMC data
  /*!
   This books histograms for the following data: AMC number, L1A number, Bunch Crossing ID, Data Length, Format Version, Run Type, Param 1, 2 and 3, Orbit Number, Board ID, GEM DAV list, Buffer Status, GEM DAV count, TTS state, Chamber Timeout, OOS GLIB sync, CRC, L1AT, and DlengthT
   */
  void bookHistograms()
  {
    m_gebsH = new GEB_histogram*[12];
    for (unsigned int i = 0; i<12; i++){
      m_gebsH[i] = 0;
    }
    m_dir->cd();
    AMCnum       = new TH1F("AMCnum", "AMC number", 12,  0, 12);
    //L1A        = new TH1F("L1A", "L1A ID", 0xffffff,  0x0, 0xffffff);
    L1A_mismatch = new TH1F("L1A_mismatch", "L1A_mismatch (AMC L1A - AMC13 LV1_idT)", 40,  -20.0, 20.0);
    BX           = new TH1F("BX", "BX ID", 4095,  0x0, 0xfff);
    //Dlength    = new TH1F("Dlength", "Data Length", 0xfffff,  0x0, 0xfffff);
    FV           = new TH1F("FV", "Format Version", 15,  0x0, 0xf);
    Rtype        = new TH1F("Rtype", "Run Type", 15,  0x0, 0xf);
    Param1       = new TH1F("Param1", "Run Param 1", 256,  -0.5, 255.5);
    Param2       = new TH1F("Param2", "Run Param 2", 256,  -0.5, 255.5);
    Param3       = new TH1F("Param3", "Run Param 3", 1024,  -0.5, 1023.5);
    mspl         = new TH1F("mspl", "mspl", 8,  -0.5, 7.5);
    Onum         = new TH1F("Onum", "Orbit Number", 0xffff,  0, 0xffff);
    BID          = new TH1F("BID", "Board ID", 0xffff,  0, 0xffff);
    GEMDAV       = new TH1F("GEMDAV", "GEM DAV list", 24,  0, 24);
    Bstatus      = new TH1F("Bstatus", "Buffer Status", 24,  0, 24);
    GDcount      = new TH1F("GDcount", "GEM DAV count", 32,  0, 32);
    Tstate       = new TH1F("Tstate", "TTS state", 15,  0, 15);
    ChamT        = new TH1F("ChamT", "Chamber Timeout", 24, 0, 24);
    OOSG         = new TH1F("OOSG", "OOS GLIB", 1, 0, 1);
  }
  
  //!Fills the histograms for AMC data
  /*!
   This fills the histograms for the following data: AMC number, L1A number, Bunch Crossing ID, Data Length, Format Version, Run Type, Param 1, 2 and 3, Orbit Number, Board ID, GEM DAV list, Buffer Status, GEM DAV count, TTS state, Chamber Timeout, OOS GLIB sync, CRC, L1AT, and DlengthT
   */
  void fillHistograms(AMCdata *amc, AMC13Event *amc13){
    AMCnum->Fill(amc->AMCnum());
    //L1A->Fill(amc->L1A());
    L1A_mismatch->Fill(amc->L1A()-amc13->LV1_idT());
    BX->Fill(amc->BX());
    //Dlength->Fill(amc->Dlength());
    FV->Fill(amc->FV());
    Rtype->Fill(amc->Rtype());
    Param1->Fill(amc->Param1());
    Param2->Fill(amc->Param2());
    Param3->Fill(amc->Param3());
    Onum->Fill(amc->Onum());
    BID->Fill(amc->BID());
    GDcount->Fill(amc->GDcount());
    Tstate->Fill(amc->Tstate());
    OOSG->Fill(amc->OOSG());
    uint8_t binFired = 0;
    for (int bin = 0; bin < 24; bin++){
      binFired = ((amc->GEMDAV() >> bin) & 0x1);
      if (binFired) GEMDAV->Fill(bin);
      binFired = ((amc->Bstatus() >> bin) & 0x1);
      if (binFired) Bstatus->Fill(bin);
      binFired = ((amc->ChamT() >> bin) & 0x1);
      if (binFired) ChamT->Fill(bin);
    }
  }
  //!Adds a GEB_histogram object to m_gebsH vector
  void addGEBH(GEB_histogram* gebH, int i){m_gebsH[i]=gebH;}
  //!Returns the m_gebsH vector
  GEB_histogram* gebsH(int i){return m_gebsH[i];}
private:
  GEB_histogram **m_gebsH;              ///<A vector of GEB_histogram
  TH1F* AMCnum;                         ///<Histogram for AMC number
  TH1F* L1A;                            ///<Histogram for L1A number
  TH1F* L1A_mismatch;                   ///<Histogram for L1A vs LV1_idT mismatch
  TH1F* BX;                             ///<Histogram for Bunch Crossing ID
  TH1F* Dlength;                        ///<Histogram for Data Length
  TH1F* FV;                             ///<Histogram for Format Version
  TH1F* Rtype;                          ///<Histogram for Run Type
  TH1F* Param1;                         ///<Histogram for Parameter 1
  TH1F* Param2;                         ///<Histogram for Parameter 2
  TH1F* Param3;                         ///<Histogram for Parameter 3
  TH1F* mspl;                         ///<Histogram for Parameter 3
  TH1F* Onum;                           ///<Histogram for Orbit Number
  TH1F* BID;                            ///<Histogram for Board ID
  TH1F* GEMDAV;                         ///<Histogram for GEM DAV list
  TH1F* Bstatus;                        ///<Histogram for Buffer Status
  TH1F* GDcount;                        ///<Histogram for GEM DAV count
  TH1F* Tstate;                         ///<Histogram for TTS state
  TH1F* ChamT;                          ///<Histogram for Chamber Timeout
  TH1F* OOSG;                           ///<Histogram for OOS GLIB
};



