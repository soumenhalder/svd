
void efficiency_calculation_4(TString fname="../Input_rootfile/data_bhabha_l4.root") {
TFile *file=new TFile("histograms/final_histograms_l4_data_bhabha.root","RECREATE");

TString a="soumen_svdResolution";
TFile *fi=new TFile(fname);
TNtuple *tre = (TNtuple*)fi->Get(a);
int nentries=tre->GetEntries();

TString b="intercept";
TNtuple *tree = (TNtuple*)fi->Get(b);
int nentries1=tree->GetEntries();
 
cout<<nentries<<" "<<nentries1<<endl;

file->cd();

int count_cluster[3][2]={0};
int count_intercept[3]={0};
int count_clustersize[10][2]={0};

Float_t L,S,PXD,SVD,u_id,v_id,U,V,U_ID,V_ID;
tree->SetBranchAddress("layer",&L);
tree->SetBranchAddress("sensor",&S);
tree->SetBranchAddress("pxd",&PXD);
tree->SetBranchAddress("svd",&SVD);
tree->SetBranchAddress("u_id",&U_ID);
tree->SetBranchAddress("v_id",&V_ID);
tree->SetBranchAddress("U",&U);
tree->SetBranchAddress("V",&V);


TH2F *UV[3][1];
for(int p=0;p<3;p++){
   TString s= "sensor";
    s+=(p+1);
    s+="uv";
    UV[p][0]= new TH2F(s,s,100,-5,5,100,-7,7);
    

}

for(int i=0;i<nentries1;i++){
   tree->GetEntry(i);
   Int_t L_=(Int_t) L;
   Int_t S_=(Int_t) S;
   Int_t PXD_=(Int_t) PXD;
   Int_t SVD_=(Int_t) SVD;
   if(!(U_ID>20 && U_ID<748 && V_ID>20 && V_ID<492)) continue;
   if(L_!=4) continue;
   UV[S-1][0]->Fill(U,V); 
   if(PXD_>0 && SVD_>0){
     count_intercept[S_-1]++;
     }
  }



float error,cls_error,intercept,cls,truehit,angle;
float size,svd,cdc,pxd,layer,sensor,side,u_id,v_id;
tre->SetBranchAddress("size",&size);
tre->SetBranchAddress("angle",&angle);
tre->SetBranchAddress("PXD",&pxd);
tre->SetBranchAddress("SVD",&svd);
tre->SetBranchAddress("CDC",&cdc);
tre->SetBranchAddress("intercept_error",&error);
tre->SetBranchAddress("cluster_position_error",&cls_error);
tre->SetBranchAddress("intercept_position",&intercept);
tre->SetBranchAddress("cluster_position",&cls);
tre->SetBranchAddress("layer",&layer);
tre->SetBranchAddress("sensor",&sensor);
tre->SetBranchAddress("side",&side);
tre->SetBranchAddress("truehit_pos", &truehit);
tre->SetBranchAddress("u_id", &u_id);
tre->SetBranchAddress("v_id", &v_id);


TH1F *residual[3][2][10];
TH1F *residual_all[3][2];
TH1F *clusterposition[3][2];
TH1F *interceptposition[3][2];
TH1F *angleuv[3][2];
TH2F *residual_clspos[3][2];
TH2F *uv_id[3][2];
for(int p=0;p<3;p++){
  for(int t=0;t<2;t++){
    for(int q=0;q<10;q++){
    TString s= "41";
    s+=(p+1);
    if(t==0) s+="V:Residual Distribution ";
    if(t==1) s+="U:Residual Distribution";
    s+="for cluster size = ";
    s+=(q+1);
    TString a,b,cc;
    a=s;b=s,cc=s;
    s+="cls-intercept";
    a+="residual_vs_clspos";
    residual[p][t][q]= new TH1F(s,s,100,-0.03,0.03);
    residual[p][t][q]->GetXaxis()->SetTitle("Residual(intercept-cluster) in cm");
    residual[p][t][q]->GetXaxis()->SetTitleOffset(0.7);
    residual[p][t][q]->GetXaxis()->SetTitleSize(0.04);
    residual[p][t][q]->SetTitle(cc);

   
   }
  TString s= "41";
    s+=(p+1);

    TString a,b,c,d,e,f,g,h,y,z,aa;
    a=s;b=s;c=s;d=s;e=s;f=s,g=s,h=s;y=s;z=s,aa=s;
    if(t==0) s+="V";
    if(t==1) s+="U"; 
   s+=":Incident angle";
   angleuv[p][t]= new TH1F(s,s,60,0,180);
   angleuv[p][t]->GetXaxis()->SetTitle("angle in degree");
   angleuv[p][t]->GetXaxis()->SetTitleOffset(0.5);
   angleuv[p][t]->GetXaxis()->SetTitleSize(0.07);
  
   
   }
   a+="U:Residual_vs_ClusterPos"; 
   b+="V:Residual_vs_ClusterPos";
   c+="clusterU";
   d+="clusterV";
   e+="interceptU";
   f+="interceptV";
   g+="u_id";
   h+="v_id";
   y+="U:Residual Distribution";
   z+="V:Residual Distribution";
   residual_clspos[p][0]= new TH2F(b,b,100,-0.03,0.03,100,-7,7);
   residual_clspos[p][1]= new TH2F(a,a,100,-0.03,0.03,100,-3,3);
   clusterposition[p][1]= new TH1F(c,c,100,-3,3);
   clusterposition[p][0]= new TH1F(d,d,100,-7,7);
   interceptposition[p][1]= new TH1F(e,e,100,-3,3);
   interceptposition[p][0]= new TH1F(f,f,100,-7,7); 
   uv_id[p][0] = new TH2F(h,h,100,-7,7,512,0,512);
   uv_id[p][1] = new TH2F(g,g,100,-4,4,768,0,768);

   residual_all[p][0]= new TH1F(z,z,100,-0.03,0.03);
   residual_all[p][1]= new TH1F(y,y,100,-0.03,0.03);
   residual_all[p][0]->GetXaxis()->SetTitle("Residual(intercept-cluster) in cm");
   residual_all[p][0]->GetXaxis()->SetTitleOffset(0.7);
   residual_all[p][0]->GetXaxis()->SetTitleSize(0.04);
   residual_all[p][1]->GetXaxis()->SetTitle("Residual(intercept-cluster) in cm");
   residual_all[p][1]->GetXaxis()->SetTitleOffset(0.7);
   residual_all[p][1]->GetXaxis()->SetTitleSize(0.04);

   residual_clspos[p][1]->GetXaxis()->SetTitleSize(0.04);
   residual_clspos[p][1]->GetXaxis()->SetTitleOffset(0.7);
   residual_clspos[p][1]->GetYaxis()->SetTitleSize(0.04);
   residual_clspos[p][1]->GetYaxis()->SetTitleOffset(0.5);
   residual_clspos[p][1]->GetXaxis()->SetTitle("Residual(intercept-cluster) in cm");
   residual_clspos[p][1]->GetYaxis()->SetTitle("Cluster position in cm");

   residual_clspos[p][0]->GetXaxis()->SetTitleSize(0.04);
   residual_clspos[p][0]->GetXaxis()->SetTitleOffset(0.7);
   residual_clspos[p][0]->GetYaxis()->SetTitleSize(0.04);
   residual_clspos[p][0]->GetYaxis()->SetTitleOffset(0.5);
   residual_clspos[p][0]->GetXaxis()->SetTitle("Residual(intercept-cluster) in cm");
   residual_clspos[p][0]->GetYaxis()->SetTitle("Cluster position in cm");
}



for(int i=0;i<nentries;i++){
   tre->GetEntry(i);
   if(!(u_id>20 && u_id<748 && v_id>20 && v_id<492)) continue;
   Int_t size_=(Int_t) size;
   Int_t sensor_ =(Int_t) sensor;
   Int_t side_ = (Int_t) side;
   if(layer!=4) continue;
   if(size_>9) count_clustersize[9][side_]++;
   else count_clustersize[size_-1][side_]++;
   if(size_>9) size_=10;
   if(!(svd>0 && pxd>0)) continue;
  
   if(TMath::Abs(cls-intercept)<0.05) count_cluster[sensor_-1][side]++;
   residual[sensor_-1][side_][size_-1]->Fill(intercept-cls);
   residual_all[sensor_-1][side_]->Fill(intercept-cls);
   angleuv[sensor_-1][side_]->Fill(angle);
   residual_clspos[sensor_-1][side_]->Fill(intercept-cls,cls);
   clusterposition[sensor_-1][side_]->Fill(cls);
   interceptposition[sensor_-1][side_]->Fill(intercept);
   if(side_==1) uv_id[sensor_-1][side_]->Fill(intercept,u_id);
   if(side_==0) uv_id[sensor_-1][side_]->Fill(intercept,v_id);
   

    
  }

gStyle->SetOptFit(kTRUE);

 
  for(int p=0;p<3;p++){
    for(int t=0;t<2;t++){
        TString s= "plots/41";
        s+=(p+1);
        TString x=s,y=s,z=s,zz=s;
        if(t==0) {s+="V:residualplot.eps"; x+="V:residual_vs_clusterpos.eps";y+="V:incident_angle.eps";zz+="V:clssize";}
        if(t==1) {s+="U:residualplot.eps"; x+="U:residual_vs_clusterpos.eps";y+="U:incident_angle.eps";zz+="U:clssize";}
        TCanvas *c1=new TCanvas("c1","c1",1000,1000);
        c1->cd();
        residual_all[p][t]->Draw();
        c1->SaveAs(s);
        residual_clspos[p][t]->Draw();
        c1->SaveAs(x);
        angleuv[p][t]->Draw();
        c1->SaveAs(y);        
      for(int q=0;q<4;q++){
         TString xyz=zz;
         xyz+=(q+1);
         xyz+=".eps";
         TCanvas *c2=new TCanvas("c2","c2",1000,1000);
         c2->cd();
         residual[p][t][q]->Draw();
         c2->SaveAs(xyz);
      }
    }
  }


  cout<<"intercept:41= "<<count_intercept[0]<<endl;
  cout<<"cluster:41V= "<<count_cluster[0][0]<<"  cluster:41U= "<<count_cluster[0][1]<<endl;
  cout<<"intercept:42= "<<count_intercept[1]<<endl;
  cout<<"cluster:42V= "<<count_cluster[1][0]<<"  cluster:42U= "<<count_cluster[1][1]<<endl;
  cout<<"intercept:43= "<<count_intercept[2]<<endl;
  cout<<"cluster:43V= "<<count_cluster[2][0]<<"  cluster:43U= "<<count_cluster[2][1]<<endl;
  for(int i=0;i<10;i++) 
     cout<<i+1<<"  V="<<count_clustersize[i][0]<<"  :  U="<<count_clustersize[i][1]<<endl;
   file->Write();
}



