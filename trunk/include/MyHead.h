
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

#include "orbitStruct.h"

using namespace std;

#define EPS 1.e-12

///////////////////////////////////// Simulation variables

//const static Int_t sky_events = 1e+9;
const static Int_t sky_events = 1e+3;
const static UInt_t random_seed = 22;

const static time_t time_stamp = time(0);       //Setting timestamp for the out files

const static TString sbi_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/map-inversion-svn/trunk/SBI_data/";
const static TString acceptance_final_plot = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/dampe-gacceptance-svn/trunk/results/1516819290_acceptance_result.root";
const static TString sbi_subsample = "010";
const static string string_sbi_subsample = "010";
const static Int_t number_SBI_files = 3;       // To be precise, at the moment of sotware writing, they are 0102800000_SBI.root 0102900000_SBI.root 0103000000_SBI.root

const static string output_log = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/map-inversion-svn/trunk/logs/";
const static string output_root = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/map-inversion-svn/trunk/results/";

////////////////////////////////////////////////////////////////////////

extern string output_path_creator(const Int_t out_choose);
extern void log_file_init(ofstream &out_file);
extern void read_SBI_data(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events,bool &good,TChain* tree,TString sbi_data_path,ofstream &out_file);
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

// --------------------------- Function used to invert the direct map: from galactic to local !
extern void sky_backtrack(Float_t sat_ra[],Float_t sat_dec[],TH2D* acc);
extern void invert_map(Double_t &costheta,Double_t &phi,Double_t l,Double_t b,Float_t sat_ra[],Float_t sat_dec[],Double_t &inv_ra,Double_t &inv_dec,Double_t inv_vector[]);
extern void from_galactic_to_celestial(Double_t &ra,Double_t &dec,Double_t l,Double_t b);
extern void from_celestial_to_local(AtPolarVect vector_out,Double_t vector_in[]);
