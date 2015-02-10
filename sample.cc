{
  gSystem->Load("libRooFit.so");
  gSystem->Load("libCommonFitting.so");

  gStyle->SetOptStat(0);

  TF1 landau("landau", "TMath::Landau(x, -10., 1.)", -20., 20.);
  TF1 gauss("gauss", "gaus(0)", -20., 20.);
  gauss.SetParameters(1., 0., 2.);

  TH1F* target = new TH1F("target", "Fit Target", 40, -20., 20.);
  TH1F* temp1 = new TH1F("temp1", "Template 1", 40, -20., 20.);
  TH1F* temp2 = new TH1F("temp2", "Template 2", 40, -20., 20.);

  target->FillRandom("landau", 2000);
  target->FillRandom("gauss", 1000);

  temp1->FillRandom("landau", 5000);

  temp2->FillRandom("gauss", 10000);

  TemplateChi2Fitter* fitter = TemplateChi2Fitter::singleton();

  fitter->setTarget(target);
  fitter->addTemplate(temp1, "temp1");
  fitter->addTemplate(temp2, "temp2");

  fitter->fit();

  double s1 = fitter->getFloat("temp1")->getVal();
  double s2 = fitter->getFloat("temp2")->getVal();

  temp1->Scale(s1);
  temp2->Scale(s2);

  temp1->SetLineColor(kRed);
  temp1->SetLineStyle(kDashed);

  temp2->SetLineColor(kBlue);
  temp2->SetLineStyle(kDashed);

  target->SetLineColor(kBlack);

  TCanvas* c1 = new TCanvas("c1", "c1");

  target->Draw();
  temp1->Draw("SAME");
  temp2->Draw("SAME");

  c1->Print("test.pdf");
}
