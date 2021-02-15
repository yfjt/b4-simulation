void angle_distribution(){
  
  TString filename = "angle_estimate.dat";
  TCanvas* c1 = new TCanvas();
  c1->Divide(5,2);
  //TH1D* hist0 = new TH1D("hist0","",120,0,3);
  TH1D* hist0 = new TH1D("hist0","angle distribution;degree;count",90,0,90);
  TCanvas* c2 = new TCanvas();
  
  
  int bin = 150;
  double d = 42; //mm
  int canvasn = 0;

  double x,y;
  double resx[3],resy[3];
  int i = 0;
  ifstream ifs0("position_resolution.dat");
  while(ifs0 >> x >> y){
    resx[i] = x;
    resy[i] = y;
    // cout << i << " " << resx[i] << " " << resy[i] << endl;
    i++;
  }
  ifs0.close();
  
  double x1,x2,x3,y1,y2,y3;
  ifstream ifs(filename);
  while(ifs >> x1 >> x2 >> x3 >> y1 >> y2 >> y3){
    c1->cd(canvasn+1);
    TGraphErrors* gx = new TGraphErrors();
    gx->SetPoint(0,x1,2*d);
    gx->SetPointError(0,resx[0],0);
    gx->SetPoint(1,x2,d);
    gx->SetPointError(1,resx[1],0);
    gx->SetPoint(2,x3,0);
    gx->SetPointError(2,resx[2],0);
    gx->SetMarkerStyle(20);
    gx->SetMarkerSize(1.0);
    gx->Draw("AP");

    double slopex;
    if(x1 != x3){
      TF1* fx = new TF1("fx","pol1",0,150);
      fx->SetParameter(0,d*x1/(x3-x1));
      fx->SetParameter(1,d/(x1-x3));
      gx->Fit("fx");
      slopex = 1/(fx->GetParameter(1));

    }else{
      slopex = 0;
    }
    
    c1->cd(canvasn+6);
    TGraphErrors* gy = new TGraphErrors();
    gy->SetPoint(0,y1,2*d);
    gy->SetPointError(0,resy[0],0);
    gy->SetPoint(1,y2,d);
    gy->SetPointError(1,resy[1],0);
    gy->SetPoint(2,y3,0);
    gy->SetPointError(2,resy[2],0);
    gy->SetMarkerStyle(20);
    gy->SetMarkerSize(1.0);
    gy->Draw("AP");
    
    double slopey;
    if(y1 != y3){
      TF1* fy = new TF1("fy","pol1",0,150);
      fy->SetParameter(0,d*y1/(y3-y1));
      fy->SetParameter(1,d/(y1-y3));
      gy->Fit("fy");
      slopey = 1/(fy->GetParameter(1));
    }else{
      slopey = 0;
    }
    
    
    double slope = sqrt(pow(slopex,2)+pow(slopey,2));
    double theta = atan(slope)*180./TMath::Pi(); //degree
    hist0->Fill(theta);
    
    //cout << slopex << " " << slopey << endl;
    canvasn += 1;
  }

  c2->cd();
  hist0->Draw("E");
  
  ifs.close();
  
  TH1D* hist1 = new TH1D("hist1","angle distribution after corection;degree;count",90,0,90);
  TCanvas* c3 = new TCanvas();
  TH1D* hist2 = new TH1D("hist2","angle distribution after corection and per solid angle;degree;count",90,0,90);
  TCanvas* c4 = new TCanvas();
  
  int j;
  double val;
  double ratio[90],in[90];
  ifstream ifs1("correction_135mm.dat");
  while(ifs1 >> j >> val){
    //k = j-1;
    ratio[j-1] = val;
    in[j-1] = hist0->GetBinContent(j);
    //cout << k << " " << ratio[k] << endl;
    hist1->SetBinContent(j,ratio[j-1]*in[j-1]);
    //cout << k << " " << ratio[k]*in[k] << endl;
  }
  c3->cd();
  hist1->Draw("E");

  TF1* fit1 = new TF1("fit1","[0]*sin(x*TMath::Pi()/180.)*pow(cos(x*TMath::Pi()/180.),[1])",0,50);
  fit1->SetParameters(100,2);
  hist1->Fit("fit1","","",0,45);
  
  ifs1.close();

  double distr[90];
  for(int k=0;k<90;k++){
    distr[k] = ratio[k]*in[k]/(2*TMath::Pi()*sin((k+0.5)*TMath::Pi()/180.));
    hist2->SetBinContent(k+1,distr[k]);
  }
  c4->cd();
  hist2->Draw("E");
  
  TF1* fit2 = new TF1("fit2","[0]*pow(cos(x*TMath::Pi()/180.),[1])",0,50);
  fit2->SetParameters(100,2);
  hist2->Fit("fit2","","",1,53);

  gStyle->SetOptFit(1111);
  
  
}
  
