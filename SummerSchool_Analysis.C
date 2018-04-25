

TCanvas * c1 = new TCanvas("c1", "Linear Scale", 800, 600);
TCanvas * c2 = new TCanvas("c2", "Log Scale", 800, 600);

// Definition of energy histogram
TH1D * h_Energy2 = new TH1D("h_Energy2", "h_Energy2", 4094, 0, 4094);
TH1D * h_EnergySpectrum = new TH1D("h_EnergySpectrum", "h_EnergySpectrum", 4094, 0, -0.210+(4094*0.532)+(0.0000000203*pow(4094,2)));


void NAA_Analysis(int GroupNum, string type){

  h_EnergySpectrum->GetXaxis()->SetTitle("Energy (keV)");
  h_EnergySpectrum->GetXaxis()->CenterTitle();
  
  // Insert name of text file that holds columns representing energy and counts
  char NAAfilename[100];
  sprintf(NAAfilename, "Spectra/SpectraPractice_%d_%s.TKA", GroupNum, type.c_str());

  int binNum = 0;
  double counts = 0.;
  double energy = 0.;

  // Open text file
  ifstream ff1;
  ff1.open(NAAfilename);

  while(!ff1.eof()){
    binNum++;
    if (binNum < 3) continue;
    ff1 >> counts ;


    energy = -0.210+(binNum*0.532)+(0.0000000203*pow(binNum,2));
    //h_EnergySpectrum->Fill(energy*0.996988 + 2.71258, counts);  // energy spectrum needed a slight correction
    h_EnergySpectrum->Fill(energy, counts);
  }
  ff1.close();


  // Define canvases for drawing energy spectrum
  //c1 = new TCanvas("c1", "Linear Scale", 800, 600);
  //c2 = new TCanvas("c2", "Log Scale", 800, 600);

  c1->cd();
  h_EnergySpectrum->Draw();
  c2->cd();
  c2->SetLogy();
  h_EnergySpectrum->Draw();
}

void FitPeak(double low, double center, double high){


  if (center < low || high < low || high < center){ 
    std::cout << "Fit parameters not in order, try again." << std::endl;
    return;
  }

  TF1 * gausfit = new TF1("gausfit", "[0]*exp(-0.5* ( (x-[1]) / [2]  )**2) + [3] +[4]*TMath::Erfc((x-[1])/[2]) ", low, high);  

  gausfit->SetParameter(0, h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(center))*2); 
  gausfit->SetParLimits(0, h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(center)), h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(center))*10);
  gausfit->SetParameter(1, center);
  gausfit->SetParLimits(1, center - 7, center + 7);
  gausfit->SetParameter(2, 5);
  gausfit->SetParLimits(2, 0, 20);
  gausfit->SetParameter(3, h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(low)));
  gausfit->SetParLimits(3, h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(low))/2, h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(low))*3);
  gausfit->SetParameter(4, h_EnergySpectrum->GetBinContent(h_EnergySpectrum->FindBin(center))/2);
  gausfit->SetParLimits(4, 0, h_EnergySpectrum->FindBin(center)*5);
  
  h_EnergySpectrum->Fit(gausfit, "r");

  TF1 * gausFit = h_EnergySpectrum->GetFunction("gausfit");
  float p0 = gausFit->GetParameter(0);
  float p1 = gausFit->GetParameter(1);
  float p2 = gausFit->GetParameter(2);
  float p3 = gausFit->GetParameter(3);
  float gausInt = p0*p2*TMath::Sqrt(6.2831853);

  cout << "Peak Position = " << p1 << endl;
  cout << "Gaussian integral = " << gausInt << endl;

}
