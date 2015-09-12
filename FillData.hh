#ifndef FILLDATA_HH
#define FILLDATA_HH

#include <string>
#include <TFile.h>
#include <TTree.h>

class FillData {
  public:
    FillData(std::string fileName="output.root");
	~FillData() {}
    void FillTree(double time, int i, int j, int k, double cI, double cB, double cR, int cn, double temperature, double quench);
    void Close();
    void Plot(std::string fName="output.root");
    double GetTheta(double xc, double yc);

  private:
    std::string fileName;
	TFile* file;
	TTree* tree;
	int* id;
	int  n;
	double I, B, R;
	double t;
	double T;
    double qchT;
};

#endif
