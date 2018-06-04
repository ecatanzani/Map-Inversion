
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <cmath>

///// ROOT libraries 

#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TFile.h"
#include "TChain.h"
#include "TColor.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TAxis.h"

#include "orbitStruct.h"

#define EPS 1.e-12
#define err_inv 1.e-3

///////////////////////////////////// Simulation variables

const static Int_t sky_events = 1e+3;
const static UInt_t random_seed = 22;

const static time_t time_stamp = time(0);       //Setting timestamp for the out files

const static TString sbi_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/Stuff/Map-Inversion/SBI_data/";
const static TString acceptance_final_plot = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/Stuff/DAMPE-GAcceptance/results/1527159147_acceptance_result.root";
const static TString sbi_subsample = "010";
const static std::string string_sbi_subsample = "010";
const static Int_t number_SBI_files = 3;       // To be precise, at the moment of sotware writing, they are 0102800000_SBI.root 0102900000_SBI.root 0103000000_SBI.root

const static std::string output_log = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/Stuff/Map-Inversion/logs/";
const static std::string output_root = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/Stuff/Map-Inversion/results/";

const static bool rndm_back_procedure = true;

////////////////////////////////////////////////////////////////////////

extern std::string output_path_creator(const Int_t out_choose);
extern void log_file_init(std::ofstream &out_file);
extern void read_SBI_data(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events,bool &good,TChain &tree,TString sbi_data_path,std::ofstream &out_file);
extern Bool_t check_sbi_loading(Float_t galactic_lat,Float_t galactic_lon,Float_t geographic_lat,Float_t geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t sec,UShort_t n_events);
extern Bool_t chech_if_null_variable(Float_t in_variable);
extern void reinitialize_all_variables(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out,Double_t right_vector[]);
extern void AtVect_To_AtPolarVect(double in_vector[],AtPolarVect &vector_out);
extern void invert_AtPolarVect_direction(AtPolarVect vector_out,AtPolarVect &vector_out_inv);
extern void AtPolarVect_to_vector(AtPolarVect &input_polar,double out_array[]);
extern void from_celestial_to_galactic(Double_t ra,Double_t dec,Double_t &l,Double_t &b);
extern void from_local_to_galactic(Double_t costheta,Double_t phi,Double_t &l,Double_t &b,Float_t sat_ra[],Float_t sat_dec[],Double_t &right_ra,Double_t &right_dec,Double_t right_vector[]);
extern void R_from_local_to_galactic(Double_t costheta,Double_t phi,Double_t &l,Double_t &b,Float_t sat_ra[],Float_t sat_dec[]);
extern void R_from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out);

// --------------------------- Function used to invert the direct map: from galactic to local !

extern void sky_backtrack(Float_t sat_ra[],Float_t sat_dec[],TH2D* evDist,TH2D *evDist_border,TH2D &direct,TH2D &inverse,TH1D &h_cos,TH1D &h_phi,TH1D &h_cos_LE,TH1D &h_phi_LE,Double_t &max_costheta,Double_t &min_costheta,Double_t &max_phi,Double_t &min_phi,Double_t tree_idx);
extern void sky_backtrack_rndm(Float_t sat_ra[],Float_t sat_dec[],TH2D *evDist_border,TH2D &inverse,Bool_t &first_call,Int_t &n_peacks,Int_t &n_valleys,std::ofstream &log_file,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b);
extern bool POI_discrete_plotting(Float_t sat_ra[],Float_t sat_dec[],TH2D* evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Bool_t &first_call,Int_t &n_peacks,Int_t &n_valleys);
extern void GetPOI(std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Float_t sat_ra[],Float_t sat_dec[],Bool_t &good_PR,Int_t n_peacks,Int_t n_valleys);
extern Double_t obtain_min_histo(TH1D *histo);
extern bool points_are_inside(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gvalley_l,Int_t n_peacks,Int_t n_valleys);
extern void obtain_galactic_max_min(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Int_t n_peacks,Int_t n_valleys,bool points_inside,Double_t &max_b,Double_t &min_b,Double_t &max_l,Double_t &min_l,Double_t &Nmax_b,Double_t &Nmin_b,Double_t &Pmax_b,Double_t &Pmin_b);
extern void R_invert_map(Double_t &costheta,Double_t &phi,Double_t l,Double_t b,Float_t sat_ra[],Float_t sat_dec[]);
extern Bool_t outside_evDist(Double_t costheta,Double_t phi,TH2D* evDist_border);
extern void get_relevant_evDist_points(TH2D* evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_phi);
extern Int_t get_bunch_dimension(TH2D* evDist_border);
extern void invert_map(Double_t &costheta,Double_t &phi,Double_t l,Double_t b,Float_t sat_ra[],Float_t sat_dec[],Double_t &inv_ra,Double_t &inv_dec,Double_t inv_vector[]);
extern void from_galactic_to_celestial(Double_t &ra,Double_t &dec,Double_t l,Double_t b);
extern void from_celestial_to_local(AtPolarVect vector_out,Double_t vector_in[]);
extern void obtain_costheta_phi(Double_t &costheta,Double_t &phi,Float_t sat_ra[],Float_t sat_dec[],Double_t vector_in[]);
extern void obtain_costheta_phi_ROOTf(Double_t &costheta,Double_t &phi,Float_t sat_ra[],Float_t sat_dec[],Double_t vector_in[]);
extern void get_acceptance_border(TH2D *histo,TH2D* histo_border);
extern Bool_t outside_acceptance(Double_t costheta,Double_t phi,TH2D *evDist_border);

