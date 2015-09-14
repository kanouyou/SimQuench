#ifndef QAnalysis_HH
#define QAnalysis_HH

#include <string>
#include <TFile.h>
#include <TTree.h>

class QAnalysis {
  public:
    QAnalysis(std::string fileName);
	~QAnalysis() {}
    void PlotTemperatureDistribution(double qtime, int numPhi, bool scale=false);
    void PlotSpotTemperature(int z, int phi, int r);
    void PlotCurrent();
	void PlotQuenchedNumber();
	void PlotResistance();
	void PlotVoltage();
	void PlotPower();
	void PlotVelocityDistribution(int phi);
	TGraph* GetTemperatureGraph(int z, int phi, int r);
    double GetMinimum(std::string name);
	double GetMaximum(std::string name);
    double GetMinTemperature(double fTime);
	double GetMaxTemperature(double fTime);
    void   PlotDeltaTemperature();

  private:
    TTree* ReadTree(TFile* file);
    
  private:
    TFile* file;
    int*   posID;
    double time;
	double current;
	double field;
	double R;
	int    qchCell;
	double temp;
	double qchTime;
};

#endif
