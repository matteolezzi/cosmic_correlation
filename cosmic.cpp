#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <math.h>
#include <time.h>
#include <unistd.h>

// ROOT libraries
#include <TMath.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TRandom3.h>

using namespace std;

// Function to calculate angular distance
// The input angles are given in radians
// Returns the angular distance in radians
double ang_dist(double ra1, double dec1, double ra2, double dec2) {
    double d;
    double PI = acos(-1.);
    d = acos(cos(PI / 2 - dec1) * cos(PI / 2 - dec2) + sin(PI / 2 - dec1) * sin(PI / 2 - dec2) * cos(ra1 - ra2));
    return d;
}

void cosmic_correlation() {
    gROOT->Reset(); // Reset ROOT
    
    int nsel = 69; // Number of selected events
    double disp = 2346.; // Scaling factor
    ifstream in;
    
    // Open the input file containing selected events
    in.open("C:/Users/Matteo Lezzi/Desktop/Matteo/Universita/Magistrale/laboratorio_di_analisi_dati/correlazione/selected69PP.txt");
    
    TTree *dataT = new TTree("dataT", "Correlation Data");
    
    // Variables for data processing
    string augerId;
    double theta = 0.; 
    double phi = 0.; 
    double Vdecl[69] = {0.};
    double Vra[69] = {0.};
    double decl = 0.;
    double ra = 0.;  
    double bgal = 0.;
    double lgal = 0.;
    double utc = 0;
    double tmp = 0;
    double eneEeV = 0.;
    double year = 0;
    int counter = 0;

    // Read data from file
    while (counter < 69 && in >> augerId >> theta >> tmp >> phi >> tmp >> lgal >> bgal >> Vra[counter] >> Vdecl[counter] >> utc >> tmp >> eneEeV) {
        decl = Vdecl[counter];
        ra = Vra[counter];  
        
        string cyear(augerId, 0, 4);
        if (!in.good()) break;
        
        cout << "List:" << "\t" << counter << "\t" << Vra[counter] << "\t" << Vdecl[counter] << "\t" << eneEeV << "\t" << endl; 
        counter++;
    }

    // Data processing and histogram creation
    double degTOrad = acos(-1.) / 180;
    double radTOdeg = 1. / degTOrad;
    TH1F *data = new TH1F("data", "data", 180, 0, 180);
    
    // Compute angular distances and fill histogram
    for (int i = 0; i < nsel; i++) {
        for (int j = i + 1; j < nsel; j++) {
            data->Fill(radTOdeg * ang_dist(degTOrad * Vra[i], degTOrad * Vdecl[i], degTOrad * Vra[j], degTOrad * Vdecl[j]));
        }
    }

    TCanvas *angc = new TCanvas("angc", "Angular Correlation", 1000, 800);
    data->Draw();
    
    // Simulation Data Processing
    int nsel_sim;
    ifstream in1;
    in1.open("C:/Users/Matteo Lezzi/Desktop/Matteo/Universita/Magistrale/laboratorio_di_analisi_dati/correlazione/eventisim69.txt");
    nsel_sim = 69000;
    
    double MC_Vdecl[69000] = {0.};
    double MC_Vra[69000] = {0.};
    
    int counter1 = 0;
    while (counter1 < 69000 && in1 >> MC_Vra[counter1] >> MC_Vdecl[counter1]) {
        counter1++;
    }
    
    // Create histogram for simulated data
    TH1F *simul = new TH1F("simul", "Simulation", 180, 0, 180);
    for (int i = 0; i < counter1 / 100; i++) {
        for (int j = i + 1; j < counter1 / 100; j++) {
            simul->Fill(radTOdeg * ang_dist(degTOrad * MC_Vra[i], degTOrad * MC_Vdecl[i], degTOrad * MC_Vra[j], degTOrad * MC_Vdecl[j]));
        }
    }
    
    simul->SetLineColor(2);
    simul->Scale(data->GetEntries() / simul->GetEntries());
    simul->Draw("same");
    
    // Significance Plot Implementation
    TCanvas *signifC = new TCanvas("significance", "Significance", 1000, 600);
    TH1F *signif = new TH1F("signifc", "Significance", 180, 0, 180);
    
    for (int i = 1; i <= 180; i++) {
        double N_data = data->GetBinContent(i);
        double N_sim = simul->GetBinContent(i);
        double significance = (N_sim > 0) ? (N_data - N_sim) / sqrt(N_sim) : 0;
        signif->SetBinContent(i, significance);
    }
    
    signif->SetLineColor(4);
    signif->Draw();
    
    // Sky Map Representation
    TCanvas *maps = new TCanvas("maps", "Sky Map", 1000, 600);
    TH2F *has = new TH2F("has", "Aitoff Projection", 360, -180, 180, 179, -89.5, 89.5);
    
    double MC_c_asc = 0;
    for (int i = 0; i < counter1; i++) {
        MC_c_asc = (MC_Vra[i] > 180) ? MC_Vra[i] - 360 : MC_Vra[i];
        has->Fill(MC_c_asc, MC_Vdecl[i]);
    }
    has->Draw("aitoff");
}
