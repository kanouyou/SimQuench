#include <iostream>
#include <string>
#include <TApplication.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "QAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
  
  TApplication* app = new TApplication("app", &argc, argv);
  
  TCanvas* c0 = new TCanvas("temp", "temp", 1000, 1000);
  c0->Divide(3,3);

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(99);

  double time = 10.;

  QAnalysis* ana = new QAnalysis("output.root");

  for (int i=0; i<8; i++) {
    c0->cd(i+1);
    gPad->SetTicks(1,1);
	gPad->SetRightMargin(0.13);
	ana->PlotTemperatureDistribution(time*(i+1), 1, true);
  }
  
  // Plot time-temperature
  TCanvas* c1 = new TCanvas("timetemp", "timetemp", 1000, 1000);
  c1->Divide(2,2);
  
  TMultiGraph* mg = new TMultiGraph();
  TGraph* gr[3];
  
  gr[0] = ana->GetTemperatureGraph( 1, 1, 18);
  gr[1] = ana->GetTemperatureGraph(58, 1, 18);
  gr[2] = ana->GetTemperatureGraph( 1, 1,  2);
  
  for (int i=0; i<3; i++) mg->Add(gr[i], "l");
  
  c1->cd(1);
  gPad->SetTicks(1,1);
  mg->SetTitle(";Time [sec];Temperature [K]");
  mg->Draw("a");
  
  // Plot Quenched Cell
  c1->cd(2);
  gPad->SetTicks(1,1);
  
  ana->PlotResistance();
  
  // plot voltage
  c1->cd(3);
  gPad->SetTicks(1,1);
  ana->PlotVoltage();
  
  // plot queched time
  c1->cd(4);
  gPad->SetTicks(1,1);
  ana->PlotVelocityDistribution(1);
  
  // Plot delta temperature
  TCanvas* c2 = new TCanvas("delta", "delta", 600, 500);
  c2->SetTicks(1,1);
  //c2->SetAttLinePS(kBlack, 1, 2);

  ana->PlotDeltaTemperature();

  app->Run();

  return 0;
}

