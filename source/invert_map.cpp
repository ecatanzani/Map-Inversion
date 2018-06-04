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
    std::string log_path = output_path_creator(0),root_out_path = output_path_creator(1);
    TString h_name;

    std::vector<Double_t> peacks_costheta,peacks_phi,valley_costheta,valley_phi;
    std::vector<Double_t> Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b;
    
    Double_t max_costheta_diff=0,min_costheta_diff=0,max_phi_diff=0,min_phi_diff=0;
    Bool_t first_call=true;
    Int_t n_peacks=0,n_valleys=0;
    
    /////////////////////////////////////////////////////////////
    
    std::ofstream output_log_file(log_path);     //log file creation !
    if(!output_log_file.is_open()) {
        std::cout<<"\n\nCannot create output file! Program finished !"<<std::endl;
        exit(-2);
    }

    log_file_init(output_log_file);

    TChain tree("SBItree");      //Defining a TTree to read SBI data fil
    read_SBI_data(g_lat,g_lon,geo_lat,geo_lon,sat_ra,sat_dec,sec,n_ev,good,tree,sbi_path,output_log_file);    //That function fills the tree reading from the SBI data files

    gRandom->SetSeed(random_seed);

    /////////////////////////////// Opening DAMPE acceptance 2D-histo
  
    TFile acc_file(acceptance_final_plot.Data());
    if(acc_file.IsZombie()) {
        std::cout << "\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
        exit(-2);
    }

    static TH2D* evDist = (TH2D*)acc_file.Get("EventsDistribution");
    
    /////////// Border of the satellite acceptance
    
    static TH2D* evDist_border = (TH2D*)evDist->Clone("evDist_border");
    evDist_border->Reset();
    get_acceptance_border(evDist,evDist_border);

    /////////////////////////////// Writing resul map root file

    TFile out_maps(root_out_path.c_str(),"RECREATE");
    if(out_maps.IsZombie()) {
        std::cout<<"\n\nError writing output ROOT file\n\n";
        exit(-2);
    }
  
    evDist_border->Write();
    
    /////////////////////////// Create histos....

    TH2D direct_map("direct_map","Direct Map (infinite statistic); #cos(#theta);  #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
    TH2D inverse_map("inverse_map","Inverse Map (infinite statistic); #cos(#theta);  #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
    TH2D R_inverse_map("inverse_map","Inverse Map (infinite statistic); #cos(#theta);  #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
    TH1D h_costheta_diff("h_costheta_diff","Costheta difference; cos(#theta_{dir}) - cos(#theta_{inv}); Entries",1000,-2e-6,2e-6);
    TH1D h_phi_diff("h_phi_diff","Phi difference; #phi_{dir}-#phi_{inv}; Entries",1000,-0.04e-3,0.04e-3);
    TH1D h_costheta_diff_LE("h_costheta_diff_LE","Costheta difference; cos(#theta_{dir}) - cos(#theta_{inv}); Entries",1000000,-0.2,0.2);
    TH1D h_phi_diff_LE("h_phi_diff_LE","Phi difference; #phi_{dir}-#phi_{inv}; Entries",1000000,-2*TMath::Pi(),2*TMath::Pi());
  
    /////////////////////////////////////////////////////////////////

    chain_entries = tree.GetEntries();
  
    for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {

        tree.GetEntry(tree_idx);

        if((sec%100000)==0)
            continue;                         //there was a bug in the SBI production
        if(!good)
            continue;                         //good second

        if (((Double_t)tree_idx/chain_entries)>(perc*0.01)) {
            std::cout<<"[ "<<perc<<" % ]"<<std::endl;
            output_log_file<<"[ "<<perc<<" % ]"<<std::endl;
            perc++;
        }

        if(g_lon>180)
            g_lon-=360;
        if(geo_lon>180)
            geo_lon-=360;

        for (int i=0; i<3; i++) {
            sat_ra[i]*=TMath::DegToRad();
            sat_dec[i]*=TMath::DegToRad();
        }

        if(rndm_back_procedure)
            sky_backtrack_rndm(sat_ra,sat_dec,evDist_border,R_inverse_map,first_call,n_peacks,n_valleys,output_log_file,peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b);
        else
            sky_backtrack(sat_ra,sat_dec,evDist,evDist_border,direct_map,inverse_map,h_costheta_diff,h_phi_diff,h_costheta_diff_LE,h_phi_diff_LE,max_costheta_diff,min_costheta_diff,max_phi_diff,min_phi_diff,tree_idx);
    
    } //end loop on seconds

    std::cout<<"\n\nMax costheta difference: "<<max_costheta_diff<<std::endl;
    std::cout<<"Min costheta difference: "<<min_costheta_diff<<std::endl;
    std::cout<<"Max phi difference: "<<max_phi_diff<<std::endl;
    std::cout<<"Min phi difference: "<<min_phi_diff<<std::endl<<std::endl;

    ////////// -------------------------------------

    //h_costheta_diff->Fit("gaus");
    //h_phi_diff->Fit("gaus");
  
    ////////////////////////// Writing final objects

    out_maps.Write();

    ////////////////////////////////////////////////
  
    std::cout<<"\n\nSimulation completed !!\n\n";
    output_log_file<<"\n\nSimulation completed !!\n\n";
  
}
