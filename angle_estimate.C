void angle_estimate(){
  srand((unsigned int)time(NULL));
  double t = 3600*10;
  int N = floor(2.1*0.8*t);
  double h = 51.; //mm
  double d = 42.;
  double incount = 0.;
  double allcount = 0.;
  int n = 9;
  double after[n];
  double in[90],out[90],obs[90];
  TF1* f =new TF1("f","(sqrt(3)*3/2)*sin(x*TMath::Pi()/180.)*pow(cos(x*TMath::Pi()/180.),2.)",0,90);
  TH2D* hist0 = new TH2D("hist0","",150,-150,300,150,-150,300);
  TH1D* hist1 = new TH1D("hist1","",90,0,90);
  TH1D* hist2 = new TH1D("hist2","",90,0,90);
  TH1D* hist3 = new TH1D("hist3","",90,0,90);
  TH1D* hist4 = new TH1D("hist4","",90,0,90);
  TH1D* hist5 = new TH1D("hist5","",90,0,90);
  TH1D* hist6 = new TH1D("hist6","",n,0,90);
  TLine* line0 = new TLine(0,0,0,150);
  TLine* line1 = new TLine(0,0,150,0);
  TLine* line2 = new TLine(0,150,150,150);
  TLine* line3 = new TLine(150,0,150,150);

  ofstream ofs("angle_estimate.dat");
  //ofstream ofs("correction.dat");
    
  auto c0 =new TCanvas("c0","c0");
  c0->SetCanvasSize(500,500);
  c0->SetWindowSize(520,520);
  TCanvas* c1 =new TCanvas();
  TCanvas* c2 =new TCanvas();
  TCanvas* c3 =new TCanvas();
  TCanvas* c4 =new TCanvas();
  
  for(int i=0;i<N;i++){
    double x0 = (rand()%15000+1)/100.; //mm
    double y0 = (rand()%15000+1)/100.;
    //cout << "x0=" << x0 <<", y0=" << y0 << endl;
    
    
    double theta0 = rand()%90; //degree
    double p = f->Eval(theta0);
    double a = (rand()%1000)/1000.;

    
    double r = h*tan(theta0*TMath::Pi()/180.);
    double theta1 = (rand()%3600)*10.; //degree
    double x3 = x0+r*cos(theta1*TMath::Pi()/180.);
    double y3 = y0+r*sin(theta1*TMath::Pi()/180.);
    double x2 = x0+(r*(h+d)/h)*cos(theta1*TMath::Pi()/180.);
    double y2 = y0+(r*(h+d)/h)*sin(theta1*TMath::Pi()/180.);
    double x1 = x0+(r*(h+2*d)/h)*cos(theta1*TMath::Pi()/180.);
    double y1 = y0+(r*(h+2*d)/h)*sin(theta1*TMath::Pi()/180.);
    
    if(a<p){    
      hist0->Fill(x2,y2);
      allcount = allcount+1;
    
      if(x1>0 && x1<=150 && y1>0 && y1<=150){
	hist1->Fill(theta0);
	incount = incount+1;
	ofs << x1 << " " << x2 << " " << x3 << " " << y1 << " " << y2 << " " << y3 << endl;
      }
    }
    if(x2>0 && x2<=150 && y2>0 && y2<=150){
      hist2->Fill(theta0);
    }
    if(x2<0 || x2>=150 || y2<0 || y2>=150){
      hist3->Fill(theta0);
    }
    
  }
  
  
  c0->cd();
  hist0->Draw("colz");
  hist0->SetStats(0);
  line0->Draw("same");
  line0->SetLineColor(2);
  line1->Draw("same");
  line1->SetLineColor(2);
  line2->Draw("same");
  line2->SetLineColor(2);
  line3->Draw("same");
  line3->SetLineColor(2);

  cout << incount/allcount << endl;
  
  c1->cd();
  hist1->Draw("E");
  hist1->SetTitle("muon's angle distribution;degree;counts");

  //c2->cd();
  //hist2->Draw("E");

  for(int j=0;j<90;j++){ 
    // in[j] = hist2->GetBinContent(j+1); 
    //out[j] = hist3->GetBinContent(j+1);
    obs[j] = hist1->GetBinContent(j+1);
    //if(in[j]>0){
    //  ratio[j] = (in[j]+out[j])/in[j];
    // }
    //else{
    //  ratio[j] = 0;
    //}
    //cout << j << " " << ratio[j] << endl;
    // hist4->SetBinContent(j,ratio[j]);
    
    //ofs << j+1 << " " << ratio[j] << endl;
  }

  int k = 0;
  int angle;
  double ratio;
  ifstream ifs("correction_135mm.dat");
  while(ifs >> angle >> ratio){
    hist5->SetBinContent(k+1,obs[k]*ratio);
    k++;
  }
  
  //c2->cd();
  //hist4->Draw("E");
  //hist4->SetTitle("correction versus degree;degree;correction");
  c3->cd();
  hist5->Draw("E");
  hist5->SetTitle("muon's angle distribution after correction;degree;counts");
  
  for(int k=0;k<n;k++){
    after[k] = hist5->Integral(90*k/n+1,90*(k+1)/n);
    hist6->SetBinContent(k+1,after[k]);
  }
 
  //int in = hist1->GetBinContent(9);
  //cout << in << endl;

  c4->cd();
  hist6->Draw("E");
  

  TF1* fit = new TF1("fit","[0]*sin(x*TMath::Pi()/180.)*pow(cos(x*TMath::Pi()/180.),[1])",0,50);
  fit->SetParameter(0,100);
  fit->SetParameter(1,2);
  hist5->Fit("fit","","",0,54);
  hist6->Fit("fit","","",0,50);

  gStyle->SetOptFit(1111);

  ofs.close();
}
