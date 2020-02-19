#include "Hardware_histogram.h"
#include "TH1.h"
#include "TH2.h"
#include <Event.h>
#define NCHANNELS 128

//!A class that creates histogram for VFAT data

class VFAT_histogram: public Hardware_histogram
{
public:
  //!Calls the base constructor
  VFAT_histogram(const std::string & filename, TDirectory * dir, const std::string & hwid):Hardware_histogram(filename, dir, hwid){}//call base constructor
  //!Calls the base constructor
  VFAT_histogram(VFAT_histogram * vH):Hardware_histogram("dummy", vH->m_dir, vH->m_HWID){}//call base constructor
  ~VFAT_histogram(){}
  
  //!Books histograms for VFAT data
  /*!
   This books histograms for the following data: Difference between crc and recalculated crc, Control Bit 1010, Control Bit 1100, Control Bit 1110, Bunch Crossing Number, Event Counter, Control Flags, and Chip ID, and Fired Channels. It also creates a subdirectory for Threshold Scans and books those histograms.
   */
  void bookHistograms(){
    m_dir->cd();
    n_hits_per_event = new TH1F("n_hits_per_event", "n_hits_per_event", 129,  -0.5 , 128.5);
    BC               = new TH1F("BC", "Bunch Crossing Number", 4096,  -0.5, 4095.5);
    BC_mismatch      = new TH1F("BC_mismatch", "Bunch Crossing VFAT-AMC (BX - BC ID) mismatch", 40,  -20.0, 20.0);
    EC               = new TH1F("EC", "Event Counter", 255,  0x0 , 0xff);
    EC_mismatch      = new TH1F("EC_mismatch", "Event Counter VFAT-AMC13 (event counter - LV1_idT) mismatch", 40,  -20.0, 20.0);
    Header           = new TH1F("Header", "Header", 32,  0x0 , 0xff);
    SlotN            = new TH1F("SlotN", "Slot Number", 24,  0, 24);
    FiredChannels    = new TH1F("FiredChannels", "FiredChannels", 128, -0.5, 127.5);
    latencyScan      = new TH1F("latencyScan",   "Latency Scan", 1024,  -0.5, 1023.5);
    latencyScan2D    = new TH2F("latencyScan2D", "Latency Scan: Chan Vs Latency", 1024, -0.5, 1023.5, NCHANNELS, -0.5, NCHANNELS-0.5);
    hitMultVsLat     = new TH2F("hitMultVsLat", "Latency Scan: Hit Multiplicity Vs Latency", 1024, -0.5, 1023.5, NCHANNELS, -0.5, NCHANNELS-0.5);
    const char *warning_labels[3] = {"Flag raised", "No channels fired", "Excessive channels fired"};
    const char *error_labels[1] = {"CRC mismatch"};
    Warnings         = new TH1I("Warnings", "Warnings", 3,  0, 3);
    for (int i = 1; i<4; i++) Warnings->GetXaxis()->SetBinLabel(i, warning_labels[i-1]);
    Errors           = new TH1I("Errors", "Critical errors", 1,  0, 1);
    for (int i = 1; i<2; i++) Errors->GetXaxis()->SetBinLabel(i, error_labels[i-1]);
    
    //get maps
    m_sn = std::stoi(m_HWID);
    readMapFromFile(m_sn,m_strip_map);
  }
  
  //!Fills histograms for VFAT data
  /*!
   This fills histograms for the following data: Difference between crc and recalculated crc, Control Bit 1010, Control Bit 1100, Control Bit 1110, Bunch Crossing Number, Event Counter, Control Flags, and Chip ID, and Fired Channels
   */
  void fillHistograms(VFATdata * vfat, AMCdata *amc, AMC13Event *amc13, long long int orbitNumber){
    BC->Fill(vfat->BC());
    BC_mismatch->Fill(vfat->BC() - amc->BX());
    EC->Fill(vfat->EC());
    EC_mismatch->Fill(vfat->EC() - amc13->LV1_idT());
    Header->Fill(vfat->Header());
    m_sn = std::stoi(m_HWID);
    //crc->Fill(vfat->crc());
    // will be returned upon reimplementation of CRCCheck FIXME
    //if (vfat->CRRCcheck()) {
    //  Errors->Fill(0);
    //}
    SlotN->Fill(m_sn);
    uint16_t chan0xf = 0;
    int n_hits_fired = 0;
    for (int chan = 0; chan < 128; ++chan) {
      if (chan < 64){
        chan0xf = ((vfat->lsData() >> chan) & 0x1);
        if(chan0xf) {
          n_hits_fired++;
          FiredChannels->Fill(chan);
        }
      } else {
        chan0xf = ((vfat->msData() >> (chan-64)) & 0x1);
        if(chan0xf) {
          n_hits_fired++;
          FiredChannels->Fill(chan);
        }
      }
    }
    n_hits_per_event->Fill(n_hits_fired);
    
  }
  void fillWarnings(){
    if (FiredChannels->GetEntries() == 0) {
      Warnings->Fill(1);
    }
    // Comment out for the moment, revision is needed later FIXME
    //else if (FiredChannels->GetEntries() > 64*b1010->GetEntries()) {
    //  Warnings->Fill(2);
    //}
  }
  //!Fills the histograms for the Threshold Scans
  void fillScanHistograms(VFATdata * vfat, int runtype, int deltaV, int latency){
    bool channelFired = false;
    int n_h_fired = 0;
    for (int i = 0; i < 128; i++){
      uint16_t chan0xf = 0;
      if (i < 64){
        chan0xf = ((vfat->lsData() >> i) & 0x1);
        if(chan0xf) {
          channelFired = true;
          latencyScan2D->Fill(latency,i);
          n_h_fired++;
        }
      } else {
        chan0xf = ((vfat->msData() >> (i-64)) & 0x1);
        if(chan0xf) {
          channelFired = true;
          latencyScan2D->Fill(latency,i);
          n_h_fired++;
        }
      }
    }// end loop on channels
    if (channelFired) {
      latencyScan->Fill(latency);
      hitMultVsLat->Fill(latency,n_h_fired);
    }
  }
  
  int * getMap(){ return m_strip_map;}
  
  
  
private:
  TH1F* n_hits_per_event;
  TH1F* BC;               ///<Histogram for Bunch Crossing Number
  TH1F* BC_mismatch;      ///<Histogram for Bunch Crossing Number mismatch between VFATs and AMC
  TH1F* EC;               ///<Histogram for Event Counter
  TH1F* EC_mismatch;      ///<Histogram for Event Counter mismatch bewteen VFATs and AMC13
  TH1F* Header;           ///<Histogram for Control Flags
  TH1F* FiredChannels;    ///<Histogram for Fired Channels (uses lsData and fmData)
  TH1F* SlotN;
  TH1F* latencyScan;
  TH2F* latencyScan2D;
  TH2F* hitMultVsLat;
  TH1I* Warnings;
  TH1I* Errors;
  int m_sn;
  int m_strip_map[128];
  
  uint16_t vfatBlockWords[12];   ///<Array of uint16_t used for setVFATBlockWords
  
  //!This puts the VFAT data in an array of uint16_t to be used for the crc check
  void setVFATBlockWords(VFATdata * vfat_)
  {
    //vfatBlockWords[11] = ((0x000f & vfat_->b1010())<<12) | vfat_->BC();
    //vfatBlockWords[10] = ((0x000f & vfat_->b1100())<<12) | ((0x00ff & vfat_->EC()) <<4) | (0x000f & vfat_->Flag());
    //vfatBlockWords[9]  = ((0x000f & vfat_->b1110())<<12) | vfat_->ChipID();
    //vfatBlockWords[8]  = (0xffff000000000000 & vfat_->msData()) >> 48;
    //vfatBlockWords[7]  = (0x0000ffff00000000 & vfat_->msData()) >> 32;
    //vfatBlockWords[6]  = (0x00000000ffff0000 & vfat_->msData()) >> 16;
    //vfatBlockWords[5]  = (0x000000000000ffff & vfat_->msData());
    //vfatBlockWords[4]  = (0xffff000000000000 & vfat_->lsData()) >> 48;
    //vfatBlockWords[3]  = (0x0000ffff00000000 & vfat_->lsData()) >> 32;
    //vfatBlockWords[2]  = (0x00000000ffff0000 & vfat_->lsData()) >> 16;
    //vfatBlockWords[1]  = (0x000000000000ffff & vfat_->lsData());
  }
  
  //!Recalculates the CRC to be compared to original CRC. Difference should be 0.
  uint16_t checkCRC(uint16_t dataVFAT[11])
  {
    uint16_t crc_fin = 0xffff;
    for (int i = 11; i >= 1; i--)
    {
      crc_fin = this->crc_cal(crc_fin, dataVFAT[i]);
    }
    return(crc_fin);
  }
  //!Called by checkCRC
  uint16_t crc_cal(uint16_t crc_in, uint16_t dato)
  {
    uint16_t v = 0x0001;
    uint16_t mask = 0x0001;
    bool d=0;
    uint16_t crc_temp = crc_in;
    unsigned char datalen = 16;
    
    for (int i=0; i<datalen; i++){
      if (dato & v) d = 1;
      else d = 0;
      if ((crc_temp & mask)^d) crc_temp = crc_temp>>1 ^ 0x8408;
      else crc_temp = crc_temp>>1;
      v<<=1;
    }
    return(crc_temp);
  }
};

