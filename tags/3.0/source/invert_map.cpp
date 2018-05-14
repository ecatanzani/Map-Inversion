////////////////////////////// Software decsription
//
//
//
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////

#include "MyHead.h"

int main(int argc,char * argv[]) {
  
  ///////////// Variables for the TTrees

  Float_t g_lat = 0,g_lon = 0,geo_lat = 0,geo_lon = 0,sat_ra[3],sat_dec[3];
  UInt_t sec = 0;
  UShort_t n_ev = 0;
  bool good = true;
  Int_t chain_entries;                  // -> Variable that stores the whole number of tchain entries

  for(Int_t idx=0; idx<3; idx++) {
    sat_ra[idx] = 0;
    sat_dec[idx] = 0;
  }
    
  //////////////////////// Variable description

  /*
    
    sat_ra[3]     -> Is the array that stores the right ascension values of the satellite
    sat_dec[3]    -> Is the array that stores the celestial declination values of the satellite

    geo_lat       -> Is the geographic latitude
    geo_lon       -> Is the geographic longitude

    g_lat         -> Is the galactic latitude
    g_lon         -> Is the galactic longitude
    
                  These two are the "absolute coordinates" used to plot the final maps !!!!

   
    n_ev          -> Is the nuber of triggered events in a second
    sec           -> Is DAMPE's acquisition second number

    good          -> Is the status of the SBI

   */

  //////////////////////////////////////
  
  ///////////////////////////////////////////////////////////// 

  Double_t perc=0;                      // -> Just used to store the percentage valu
  string log_path = output_path_creator(0),root_out_path = output_path_creator(1);
  TString h_name;

  Double_t max_costheta_diff,min_costheta_diff,max_phi_diff,min_phi_diff;
  
  ofstream output_log_file(log_path);     //log file creation !
  if(!output_log_file.is_open()) {
    cout<<"\n\nCannot create output file! Program finished !"<<endl;
    exit(-1);
  }

  log_file_init(output_log_file);

  TChain *tree= new TChain("SBItree");      //Defining a TTree to read SBI data fil
  read_SBI_data(g_lat,g_lon,geo_lat,geo_lon,sat_ra,sat_dec,sec,n_ev,good,tree,sbi_path,output_log_file);    //That function fills the tree reading from the SBI data files

  gRandom->SetSeed(random_seed);

  /////////////////////////////// Opening DAMPE acceptance 2D-histo
  
TFile* acc_file= new TFile(acceptance_final_plot.Data());
  if(acc_file->IsZombie()) {
    cout << "\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
    output_log_file << "\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
    exit(-1);
  }

  static TH2D* acc = (TH2D*)acc_file->Get("Acceptance");
    
    /////////// Border of the satellite acceptance
    
    static TH2D *acc_border = (TH2D*)acc->Clone("acc_border");
    acc_border->Reset();
    get_acceptance_border(acc,acc_border);

    
  /////////////////////////////// Writing resul map root file

  TFile *out_maps = new TFile(root_out_path.c_str(),"RECREATE");
  if(out_maps->IsZombie()) {
    cout<<"\n\nError writing output ROOT file\n\n";
    exit(-1);
  }
  
        acc_border->Write();
    
    
  /////////////////////////// Create histos....

  TH2D* direct_map = new TH2D("direct_map","Direct Map (infinite statistic); #cos(#theta);  #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
  TH2D* inverse_map = new TH2D("inverse_map","Inverse Map (infinite statistic); #cos(#theta);  #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
  TH1D* h_costheta_diff = new TH1D("h_costheta_diff","Costheta difference; cos(#theta_{dir}) - cos(#theta_{inv}); Entries",1000,-2e-6,2e-6);
  TH1D* h_phi_diff = new TH1D("h_phi_diff","Phi difference; #phi_{dir}-#phi_{inv}; Entries",1000,-0.04e-3,0.04e-3);
  TH1D* h_costheta_diff_LE = new TH1D("h_costheta_diff_LE","Costheta difference; cos(#theta_{dir}) - cos(#theta_{inv}); Entries",1000000,-0.2,0.2);
  TH1D* h_phi_diff_LE = new TH1D("h_phi_diff_LE","Phi difference; #phi_{dir}-#phi_{inv}; Entries",1000000,-2*TMath::Pi(),2*TMath::Pi());
  
  /////////////////////////////////////////////////////////////////

  chain_entries = tree->GetEntries();
  
  for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {

    tree->GetEntry(tree_idx);

    if((sec%100000)==0)
      continue;                         //there was a bug in the SBI production
    if(!good)
      continue;                         //good second

    if (((Double_t)tree_idx/chain_entries)>(perc*0.01)) {
      cout<<"[ "<<perc<<" % ]"<<endl;
      output_log_file<<"[ "<<perc<<" % ]"<<endl;
      perc++;
    }

    if (g_lon>180)
      g_lon-=360;
    if (geo_lon>180)
      geo_lon-=360;

    for (int i=0; i<3; i++) {
      sat_ra[i]*=TMath::DegToRad();
      sat_dec[i]*=TMath::DegToRad();
    }

    sky_backtrack(sat_ra,sat_dec,acc,acc_border,direct_map,inverse_map,h_costheta_diff,h_phi_diff,h_costheta_diff_LE,h_phi_diff_LE,max_costheta_diff,min_costheta_diff,max_phi_diff,min_phi_diff,tree_idx);
    
  } //end loop on seconds

  cout<<"\n\nMax costheta difference: "<<max_costheta_diff<<endl;
  cout<<"Min costheta difference: "<<min_costheta_diff<<endl;
  cout<<"Max phi difference: "<<max_phi_diff<<endl;
  cout<<"Min phi difference: "<<min_phi_diff<<endl<<endl;

  ////////// -------------------------------------

  //h_costheta_diff->Fit("gaus");
  //h_phi_diff->Fit("gaus");
  
  ////////////////////////// Writing final objects

  out_maps->Write();

  ////////////////////////////////////////////////
  
  cout<<"\n\nSimulation completed !!\n\n";
  output_log_file<<"\n\nSimulation completed !!\n\n";
  
}
