struct dqmerror
{
  unsigned int errorCode;
  std::string amcID;
  std::string chamberID;
  std::string chipID;
  std::string errorMessage;
}
struct EventError
{
  long int eventNumber;
  std::vector<dqmerror> errors;
};

class GEMDQMchecker
{
  public:
    GEMDQMchecker(int eventNumber_)
    {
      //m_error = {-1,0,"N/A","N/A","N/A","N/A"};
      m_error.eventNumber = eventNumber_;
      m_error.errors.clear();
      m_eventStatus = true;
      m_checkPerformed = false;
    }
    virtual ~GEMDQMchecker(){}
    void setEventStatus(bool eventStatus_){m_eventStatus = eventStatus_;}
    bool getEventStatus()
    {
      if (m_checkPerformed) return m_eventStatus;
      else throw std::logic_error("Event status check should be completed before calling this method");
    }
    void setCheckStatus(bool checkPerformed){m_checkPerformed = checkPerformed;}
    ErrorDetails getError(){return m_error;}
  private:
    bool m_eventStatus;
    bool m_checkPerformed;
    EventError m_error;
};
