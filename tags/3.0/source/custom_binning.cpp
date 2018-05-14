
#include "MyHead.h"

void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale) {
  binning = new Double_t[n_bin_lat+1];

  if (cos_scale==false) {
    binning[0]=lat_bin_min;
    for(Int_t idx_b=1; idx_b<=n_bin_lat; idx_b++)
      binning[idx_b] = binning[idx_b-1] + (lat_bin_max-lat_bin_min)/n_bin_lat;
  }
  else {
    Double_t binning_shift = (TMath::Cos(TMath::DegToRad()*(lat_bin_max-90.0))-TMath::Cos(TMath::DegToRad()*(lat_bin_min-90.0)))/n_bin_lat;  //I subtrack 90 degrees to refer the costheta to the y axis !!
    Double_t prev_bin;  //Define a new variable for the "previous" bin value
    binning[0]=lat_bin_min;
    for(Int_t idx_b=1; idx_b<=n_bin_lat; idx_b++) {
      prev_bin = TMath::DegToRad()*(binning[idx_b-1]-90.0);
      binning[idx_b] = 90.0-(TMath::RadToDeg()*TMath::ACos(TMath::Cos(prev_bin)+binning_shift));
    }
  }
}
