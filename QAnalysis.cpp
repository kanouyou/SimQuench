#include <iostream>
#include <TH2F.h>
#include <TGraph.h>
#include "QAnalysis.h"

using namespace std;

QAnalysis::QAnalysis(string fileName) {
  file = new TFile(fileName.c_str());
}

TTree* QAnalysis::ReadTree(TFile* file) {
  TTree* tree = (TTree*)file->Get("tree");
  
  posID = new int[3];
  tree->SetBranchAddress(         "id",    posID);
  tree->SetBranchAddress(       "time",    &time);
  tree->SetBranchAddress(    "current", &current);
  tree->SetBranchAddress(      "field",   &field);
  tree->SetBranchAddress( "resistance",       &R);
  tree->SetBranchAddress( "quenchcell", &qchCell);
  tree->SetBranchAddress("temperature",    &temp);
  tree->SetBranchAddress( "quenchTime", &qchTime);

  return tree;
}

void QAnalysis::PlotTemperatureDistribution(double qtime, int numPhi, bool scale) {
  int zbin = 60;
  int rbin = 19;
  TH2F*  hist = new TH2F(Form("tempDis%.1f", qtime), Form("tempDis%.1f", qtime), zbin, 0, 60, rbin, 0, 19);
  TTree* tree = ReadTree(file);
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( time==qtime && posID[1]==numPhi ) 
      hist->Fill(posID[0], posID[2], temp);
  }
  
  if (scale==true) hist->GetZaxis()->SetRangeUser(0., GetMaximum("temp"));
  hist->SetTitle(Form("time = %.1f [sec]; Z; R; Temperature [K]", qtime));
  hist->Draw("colz");
}

void QAnalysis::PlotSpotTemperature(int z, int phi, int r) {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==z && posID[1]==phi && posID[2]==r ) {
	  gr->SetPoint(cnt, time, temp);
	  cnt += 1;
    }
  }

  gr->SetTitle(Form("Position: [%i, %i, %i]; Time [sec]; Temperature [K]", z, phi, r));
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  gr->Draw("al");

  cout << "Finished PlotSpotTemperature()" << endl;
}

TGraph* QAnalysis::GetTemperatureGraph(int z, int phi, int r) {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==z && posID[1]==phi && posID[2]==r ) {
	  gr->SetPoint(cnt, time, temp);
	  cnt += 1;
    }
  }

  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  return gr;
}

void QAnalysis::PlotCurrent() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int z = 1; int phi = 1; int r = 18;
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==z && posID[1]==phi && posID[2]==r ) {
      gr->SetPoint(cnt, time, current);
	  cnt += 1;
	}
  }
  
  gr->SetTitle("Current Decay; Time [sec]; Current [A]");
  gr->SetLineColor(kAzure+7);
  gr->SetLineWidth(2);
  gr->Draw("al"); 
}

double QAnalysis::GetMaximum(string name) {
  TTree* tree = ReadTree(file);
  double max  = tree->GetMaximum(name.c_str());
  return max;
}

double QAnalysis::GetMinimum(string name) {
  TTree* tree = ReadTree(file);
  double min  = tree->GetMinimum(name.c_str());
  return min;
}

void QAnalysis::PlotQuenchedNumber() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int meshNum = 60 * 19 * 4;
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==1 && posID[1]==1 && posID[2]==18 ) {
      gr->SetPoint(cnt, time, static_cast<double>(qchCell)/meshNum);
	  cnt += 1;
	}
  }

  gr->SetTitle("Quenched Cell; Time [sec]; Quenched Ratio [%]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kBlue+1);
  gr->Draw("al");
}

void QAnalysis::PlotResistance() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==1 && posID[1]==1 && posID[2]==18 ) {
      gr->SetPoint(cnt, time, R);
	  cnt += 1;
	}
  }

  gr->SetTitle("; Time [sec]; Resistance [#Omega]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kBlue+1);
  gr->Draw("al");
}

void QAnalysis::PlotVoltage() {
  TTree* tree = ReadTree(file);
  TGraph* gr  = new TGraph();
  
  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if ( posID[0]==1 && posID[1]==1 && posID[2]==18 ) {
      gr->SetPoint(cnt, time, R * current);
	  cnt += 1;
	}
  }

  gr->SetTitle("; Time [sec]; Voltage [V]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kPink-1);
  gr->Draw("al");
}

void QAnalysis::PlotPower() {
  TTree* tree = (TTree*)ReadTree(file);
  TGraph* gr  = new TGraph();

  int cnt = 0;
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
	if (posID[0]==1 && posID[1]==1 && posID[2]==18) {
      gr->SetPoint(cnt, time, R*current*current);
	  cnt++;
	}
  }
  
  gr->SetTitle(" Time [sec]; Power [W]");
  gr->SetLineWidth(2);
  gr->SetLineColor(kOrange);
  gr->Draw("al");
}

void QAnalysis::PlotVelocityDistribution(int phi) {
  TTree* tree = ReadTree(file);
  
  int zbin = 60; int rbin = 19;
  TH2F*  hist = new TH2F("quenchv", "quenchv", zbin, 0, 60, rbin, 0, 19);
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if ( time==70. && posID[1]==phi )
	  hist->Fill(posID[0], posID[2], qchTime);
  }
  
  hist->SetTitle("; Z; R; Quenched Time [sec]");
  hist->Draw("colz");
}

double QAnalysis::GetMaxTemperature(double fTime) {
  TTree* tree = (TTree*)ReadTree(file);
  vector<double> fTemperature;

  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (time==fTime) 
      fTemperature.push_back(temp);
  }
  
  double fMaximum = fTemperature[0];

  for (vector<int>::size_type i=0; i<fTemperature.size(); i++) {
    if (fMaximum<fTemperature[i])
	  fMaximum = fTemperature[i];
  }

  return fMaximum;
}

double QAnalysis::GetMinTemperature(double fTime) {
  TTree* tree = (TTree*)ReadTree(file);
  vector<double> fTemperature;
  
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (time==fTime) 
      fTemperature.push_back(temp);
  }
  
  double fMinimum = fTemperature[0];
  cout << fTemperature.size() << endl;

  for (vector<int>::size_type i=0; i<fTemperature.size(); i++) {
    if (fMinimum>fTemperature[i])
	  fMinimum = fTemperature[i];
  }

  return fMinimum;
}

void QAnalysis::PlotDeltaTemperature() {
  TGraph* gr = new TGraph();
  int nt = 80;
  double fMin, fMax;
  
  for (int i=0; i<nt; i++) {
    fMin = GetMinTemperature(static_cast<double>(i+1));
	fMax = GetMaxTemperature(static_cast<double>(i+1));
    gr->SetPoint(i, static_cast<double>(i+1), fMax-fMin);
	cout << "Time: " << i << " [sec]; Maximum: " << fMax << " [K]; Minimum: " << fMin << " [K]" << endl;
  }

  gr->SetLineStyle(2);
  gr->SetLineWidth(2);
  gr->SetLineColor(kPink-2);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kPink-2);
  gr->SetTitle("; Time [sec]; #Delta T [K]");
  gr->Draw("apl");
}
