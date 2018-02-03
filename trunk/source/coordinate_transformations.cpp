
#include "MyHead.h"

void from_local_to_galactic(Double_t costheta,Double_t phi,Double_t &l,Double_t &b,Float_t sat_ra[],Float_t sat_dec[],Double_t &right_ra,Double_t &right_dec,Double_t right_vector[]) {
  Double_t ra=-999, dec=-999;
  AtVect vector_in;
  AtPolarVect vector_out, vector_out_inv;

  vector_in[0]=sin(acos(costheta))*cos(phi);
  vector_in[1]=sin(acos(costheta))*sin(phi);
  vector_in[2]=costheta;
  from_satellite_to_celestial(sat_ra,sat_dec,vector_in,vector_out,right_vector);
  invert_AtPolarVect_direction(vector_out,vector_out_inv);
  ra=vector_out_inv.lon;
  dec=vector_out_inv.lat;
  ra=ra*TMath::RadToDeg();
  dec=dec*TMath::RadToDeg();

  right_ra=ra;
  right_dec=dec;
  
  from_celestial_to_galactic(ra,dec,l,b);
}

void from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out,Double_t right_vector[]) {
  Float_t ux1[3];
  Float_t uy1[3];
  Float_t uz1[3];
  Float_t rax = ra[0];
  Float_t ray = ra[1];
  Float_t raz = ra[2];
  Float_t decx = dec[0];
  Float_t decy = dec[1];
  Float_t decz = dec[2];

  AtVect tmp_vector_out;

  ux1[0] = cos(decx)*cos(rax);
  ux1[1] = cos(decx)*sin(rax);
  ux1[2] = sin(decx);

  uy1[0] = cos(decy)*cos(ray);
  uy1[1] = cos(decy)*sin(ray);
  uy1[2] = sin(decy);

  uz1[0] = cos(decz)*cos(raz);
  uz1[1] = cos(decz)*sin(raz);
  uz1[2] = sin(decz);
  
  for (Int_t idx=0; idx<3; idx++) {
    tmp_vector_out[idx] = vectorin[0] * ux1[idx] + vectorin[1] * uy1[idx] + vectorin[2] * uz1[idx];
    right_vector[idx]=tmp_vector_out[idx];
  }
  
  AtVect_To_AtPolarVect(tmp_vector_out,vector_out);

}

void AtVect_To_AtPolarVect(double in_vector[],AtPolarVect &vector_out) {
  double  norm01, c, s;

  norm01 = TMath::Power(in_vector[0],2) + TMath::Power(in_vector[1],2);
  
  if ( ( vector_out.r = sqrt( norm01 + TMath::Power(in_vector[2],2) ) ) == 0.0 )
    vector_out.lon = vector_out.lat = 0.0;

  norm01 = sqrt(norm01);
  vector_out.lat = asin(in_vector[2]/vector_out.r);
  c = in_vector[0]/norm01;
  s = in_vector[1]/norm01;

  if(norm01 < EPS)
    vector_out.lon = 0.0;
  else if (fabs(s) < EPS)
    vector_out.lon = (1.0 - c/fabs(c))*TMath::Pi()/2.0;
   
  else
    vector_out.lon = atan((1.0-c)/s)*2.0;

  while(vector_out.lon >= 2.*TMath::Pi())
    vector_out.lon -= 2.*TMath::Pi();

  while(vector_out.lon < 0.0)
    vector_out.lon += 2.*TMath::Pi();  

}

void invert_AtPolarVect_direction(AtPolarVect vector_out,AtPolarVect &vector_out_inv) {
  AtVect in_array;
  AtVect out_array;

  AtPolarVect_to_vector(vector_out,in_array);
  for(Int_t idx=0; idx<3; idx++)
    out_array[idx]=-in_array[idx];
  AtVect_To_AtPolarVect(out_array,vector_out_inv);
}

void AtPolarVect_to_vector(AtPolarVect &input_polar,double out_array[]) {
  out_array[0] = (input_polar.r)*cos(input_polar.lat)*cos(input_polar.lon);
  out_array[1] = (input_polar.r)*cos(input_polar.lat)*sin(input_polar.lon);
  out_array[2] = (input_polar.r)*sin(input_polar.lat);
}

void from_celestial_to_galactic(Double_t ra,Double_t dec,Double_t &l,Double_t &b) {
  Double_t ragc=192.85948,decgc=27.12825,lcp=122.932;

  Double_t ragcr=ragc*TMath::DegToRad();
  Double_t decgcr=decgc*TMath::DegToRad();
  Double_t lcpr=lcp*TMath::DegToRad();
  Double_t decr=dec*TMath::DegToRad();
  Double_t rar=ra*TMath::DegToRad();

  Double_t br=asin(sin(decgcr)*sin(decr)+cos(decgcr)*cos(decr)*(cos(rar-ragcr)));

  Double_t sin_t=(cos(decr)*sin(rar-ragcr)/cos(br));
  Double_t cos_t=((cos(decgcr)*sin(decr)-sin(decgcr)*cos(decr)*cos(rar-ragcr))/cos(br));

  Double_t t=0;
  if((sin_t>=0)&&(cos_t>=0))
    t=asin(sin_t);
  else if((sin_t>=0)&&(cos_t<0))
    t=TMath::Pi()-asin(sin_t);
  else if((sin_t<0)&&(cos_t<0))
    t=TMath::Pi()-asin(sin_t);
  else if((sin_t<0)&&(cos_t>=0))
    t=asin(sin_t);

  Double_t lr = lcpr-t;

  if(lr<0) {lr=lr+(2.*TMath::Pi());}
  if(lr>(2.*TMath::Pi())){lr=lr-(2.*TMath::Pi());}

  b = br*TMath::RadToDeg();
  l = lr*TMath::RadToDeg();

}
