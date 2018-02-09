
#include "MyHead.h"

void read_SBI_data(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events,bool &good,TChain* tree,TString sbi_data_path,ofstream &out_file) {

  // Bool_t chain_status=true;
  Bool_t chain_single_tree_check=true;
  Int_t chain_entries=(Int_t)tree->GetEntries(),bad_tree_number=-1;
  
  tree->Add(Form("%s/%s*_SBI.root", sbi_path.Data(),sbi_subsample.Data()));

  tree->SetBranchAddress("second",&sec);
  tree->SetBranchAddress("goodsbi",&good);
  tree->SetBranchAddress("glon",&galactic_lon);
  tree->SetBranchAddress("glat",&galactic_lat);
  tree->SetBranchAddress("nev",&n_events);
  tree->SetBranchAddress("ra_scx",&sat_ra[0]);
  tree->SetBranchAddress("ra_scy",&sat_ra[1]);
  tree->SetBranchAddress("ra_scz",&sat_ra[2]);
  tree->SetBranchAddress("dec_scx",&sat_dec[0]);
  tree->SetBranchAddress("dec_scy",&sat_dec[1]);
  tree->SetBranchAddress("dec_scz",&sat_dec[2]);
  tree->SetBranchAddress("lat_geo",&geographic_lat);
  tree->SetBranchAddress("lon_geo",&geographic_lon);

  tree->GetEntry(0);  //Set initial value for the variables

  /*
  
  for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {
    tree->GetEntry(tree_idx);
    chain_status *= check_sbi_loading(galactic_lat,galactic_lon,geographic_lat,geographic_lon,sat_ra,sat_dec,sec,n_events);
    if(!chain_status)
      break;
    else
      reinitialize_all_variables(galactic_lat,galactic_lon,geographic_lat,geographic_lon,sat_ra,sat_dec,sec,n_events);
  }

  */
    
  for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {
    tree->GetEntry(tree_idx);
    if(tree->GetTreeNumber()==bad_tree_number) {
      chain_single_tree_check=false;
      break;
    }
  }

  // if(chain_status==true && chain_single_tree_check==true) {
  if(chain_single_tree_check) {
    cout<<"\nSBI data has been loaded"<<endl<<endl;
    out_file<<"SBI data has been loaded"<<endl<<endl;
  }
  else {
    cout<<"\nError loading SBI data !!!"<<endl<<endl;
    out_file<<"Error loading SBI data !!!"<<endl<<endl;
    exit(-1);
  }

  /*
    ATTENCTION: just a comment on the 2 check performed

    To validate the authenticity of the TChain two different test were performed. 
    The first one initilize each Float_t variable to an initial value of 0. Then a loop on the TChain entries is performed and the function "check_sbi_loading" checks if that initial values were chainged (as they should be).
    You could say that maybe a value could correctly be 0 at a such second, for example the triggered events. Infact if a variable maybe will be "corretly" 0, the check function'll flag that event as an error in loading the SBI files.
    So that first method could not be used.

    The second one, at the countrary, is good.
    Me create a variable "bad_tree_number" with a number that could not me the real read TTree number inside the TChain.
    In that way we can check is the TCain is badly filled.

  */
}

Bool_t check_sbi_loading(Float_t galactic_lat,Float_t galactic_lon,Float_t geographic_lat,Float_t geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t sec,UShort_t n_events) {
  Bool_t status=true;
  vector<Float_t> variables;
  Int_t idx_v=0;
  
  variables.push_back(galactic_lat);
  variables.push_back(galactic_lon);
  variables.push_back(geographic_lat);
  variables.push_back(geographic_lon);
  for(Int_t idx=0; idx<3; idx++) {
    variables.push_back(sat_ra[idx]);
    variables.push_back(sat_dec[idx]);
  }
  variables.push_back(sec);
  variables.push_back(n_events);

  // for(Int_t idx_v=0; idx_v<(variables.size()); idx_v++)
  // cout<<endl<<variables.at(idx_v);

  while(status==true && idx_v<variables.size()) { 
    status*=chech_if_null_variable(variables.at(idx_v));
    idx_v++;
  }
  
  return status;
}

Bool_t chech_if_null_variable(Float_t in_variable) {
  Bool_t is_null;

  if(in_variable==0)
    is_null=true;
  else
    is_null=false;

  return !is_null;
   
}

void reinitialize_all_variables(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events) {

  galactic_lat=0;
  galactic_lon=0;
  geographic_lat=0;
  geographic_lon=0;

  for(Int_t idx=0; idx<3; idx++) {
    sat_ra[idx]=0;
    sat_dec[idx]=0;
  }

  sec=0;
  n_events=0;

}
