
#include "MyHead.h"

///////////////////////////////////////////// Math custom functions

static Double_t acot(Double_t value) { return atan(1/value); }

//////////////////////////////////////////////////////////////////////////////

void sky_backtrack(Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH2D* direct,TH2D* inverse,TH1D* h_cos,TH1D* h_phi,TH1D* h_cos_LE,TH1D* h_phi_LE,Double_t &max_costheta,Double_t &min_costheta,Double_t &max_phi,Double_t &min_phi,Double_t tree_idx) {
  
  Double_t costheta=-999,phi=-999;
  Double_t b=0,l=0;

  /////////// Variables for inverse map checkin
  
  Double_t right_ra,right_dec,inv_ra,inv_dec;
  Double_t right_vector[3],inv_vector[3];
  Double_t right_costheta,right_phi;
  Double_t diff_costheta,diff_phi;
  
  //////////////////////////////////////////////
  
  for(Int_t idx_t=0; idx_t<sky_events; idx_t++) {
    
    right_ra=-999;
    right_dec=-999;
    inv_ra=-999;
    inv_dec=-999;
    
    acc->GetRandom2(costheta,phi);

    right_costheta=costheta;
    right_phi=phi;

    from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec,right_ra,right_dec,right_vector);  //Ok, now I have the coordinates into the galactic frame. I try to go back returning to the local frame.

    //////////// Control for "nan" longitude values values
    if(isnan(l)) {
      idx_t--;
      continue;
    }

    direct->Fill(costheta,phi);
    
    invert_map(costheta,phi,l,b,sat_ra,sat_dec,inv_ra,inv_dec,inv_vector);

    inverse->Fill(costheta,phi);
    
    diff_costheta=right_costheta-costheta;
    diff_phi=right_phi-phi;
    
    h_cos->Fill(diff_costheta);
    h_phi->Fill(diff_phi);
    h_cos_LE->Fill(diff_costheta);
    h_phi_LE->Fill(diff_phi);
    
    if(tree_idx==0) {
      max_costheta=min_costheta=diff_costheta;
      max_phi=min_phi=diff_phi;
    }
    else {
      if(diff_costheta>max_costheta)
	max_costheta=diff_costheta;
      if(diff_costheta<min_costheta)
	min_costheta=diff_costheta;
      if(diff_phi>max_phi)
	max_phi=diff_phi;
      if(diff_phi<min_phi)
	min_phi=diff_phi;
    }
    
    
    ////////////////// Check tor correct map inversion

    /*
    if(fabs(right_ra-inv_ra)>err_inv)
      cout<<"\nError inverting the map !! Right ra: "<<right_ra<<" Inverted ra: "<<inv_ra<<"\t-> Diff: "<<right_ra-inv_ra<<endl;
    if(fabs(right_dec-inv_dec)>err_inv)
      cout<<"\nError inverting the map !! Right dec: "<<right_dec<<" Inverted dec: "<<inv_dec<<endl;

    for(Int_t idx=0; idx<3; idx++)
      if(fabs(right_vector[idx]-inv_vector[idx])>err_inv) {
	  cout<<"\n\nError inverting the map !"<<endl;
	  cout<<"Right "<<idx<<" component: "<<right_vector[idx]<<" inv "<<idx<<" component: "<<inv_vector[idx];
      }
  
    
    if(fabs(right_costheta-costheta)>err_inv)
      cout<<"\nError inverting the map: right costheta: "<<right_costheta<<" inv costheta: "<<costheta;
    if(fabs(right_phi-phi)>err_inv)
      cout<<"\nError inverting the map: right phi: "<<right_phi<<" inv phi: "<<phi;
    */

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    
  } //end for on "sky-events"

}

void invert_map(Double_t &costheta,Double_t &phi,Double_t l,Double_t b,Float_t sat_ra[],Float_t sat_dec[],Double_t &inv_ra,Double_t &inv_dec,Double_t inv_vector[]) {
  Double_t ra=-999, dec=-999;
  AtPolarVect vector_out, vector_out_inv;
  AtVect vector_in;
  
  from_galactic_to_celestial(ra,dec,l,b);
  //cout<<"\n\nAfter inverted: ra: "<<ra<<"\tdec: "<<dec<<endl;
  
  inv_ra=ra;
  inv_dec=dec;

  ////// Switched back to rad 
  ra=ra/TMath::RadToDeg();
  dec=dec/TMath::RadToDeg();
  ////// Build the "inverted" out struct
  vector_out_inv.lon=ra;
  vector_out_inv.lat=dec;
  vector_out_inv.r=1;
  ////// Back to the "non-inverted" out struct
  invert_AtPolarVect_direction(vector_out_inv,vector_out);
  ////// Now back to local satellite coordinates !!!
  from_celestial_to_local(vector_out,vector_in);
  
  for(Int_t idx=0; idx<3; idx++)
    inv_vector[idx]=vector_in[idx];

  ///// Now obtain costheta and phi of the icnoming CR to the satellite

  obtain_costheta_phi(costheta,phi,sat_ra,sat_dec,vector_in);

}


void from_galactic_to_celestial(Double_t &ra,Double_t &dec,Double_t l,Double_t b) {
  
  Double_t ragc=192.85948,decgc=27.12825,lcp=122.932;
  
  Double_t ragcr=ragc*TMath::DegToRad();
  Double_t decgcr=decgc*TMath::DegToRad();
  Double_t lcpr=lcp*TMath::DegToRad();
  Double_t decr,rar;
  
  Double_t old_ra=ra;
  Double_t tmp_l=0,tmp_b=0;
  
  Double_t br,lr,sin_t,cos_t,t;

  ////////////////////////////////////////////////

  br = b/TMath::RadToDeg();
  lr = l/TMath::RadToDeg();
  
  t = lcpr-lr;
  sin_t = sin(t);
  cos_t = cos(t);

  //////////// Constants to easily invert the map
  
  Double_t c1 = sin(br);
  Double_t c2 = sin(decgcr);
  Double_t c3 = cos(decgcr);
  Double_t c4 = sin_t*cos(br); 
  Double_t c5 = cos_t*cos(br); 
  Double_t c6 = ragcr;
  
  ////////////////////////////////////////////////
  
  decr = asin(c1/c2-((c3*c4)/c2)*((c3*c1-c5*c2)/(TMath::Power(c3,2)*c4+TMath::Power(c2,2)*c4)));
  rar = acot((c3*c1-c5*c2)/(TMath::Power(c3,2)*c4+TMath::Power(c2,2)*c4))+c6;
  
  ///////////////////////////////////////////////
  
  ra = rar/TMath::DegToRad();
  dec = decr/TMath::DegToRad();

  //////////////////////////////// BLACK MAGIC TO SOLVE 180 DEGREES BIAS !!!
  
  from_celestial_to_galactic(ra,dec,tmp_l,tmp_b);
  if(fabs(tmp_l-l)>err_inv || fabs(tmp_b-b)>err_inv)
    ra+=180;
  from_celestial_to_galactic(ra,dec,tmp_l,tmp_b);
  if(fabs(tmp_l-l)>err_inv || fabs(tmp_b-b)>err_inv)
    ra=old_ra-180;

  while(ra>360)
    ra-=360;

  /////////////////////////////////////////////////////////////////////////////////////////
  
}

void from_celestial_to_local(AtPolarVect vector_out,Double_t vector_in[]) {
  Double_t abs_s,s;
  Double_t c,norm01;

  norm01 = TMath::Power(vector_out.r,2)-TMath::Power(vector_out.r*sin(vector_out.lat),2);
  c=(1-TMath::Power(tan(vector_out.lon/2.),2))/(1+TMath::Power(tan(vector_out.lon/2.),2));
  abs_s = sqrt(1-TMath::Power(c,2));
  s=((1-c)/tan(vector_out.lon/2.));
  
  if(abs_s>EPS) {
    vector_in[0]=c*sqrt(norm01);
    vector_in[1]=s*sqrt(norm01);
    vector_in[2]=vector_out.r*sin(vector_out.lat);
  }
  else {
    c=1;
    vector_in[0]=c*sqrt(norm01);
    vector_in[1]=0;
    vector_in[2]=vector_out.r*sin(vector_out.lat);
  }

}

void obtain_costheta_phi(Double_t &costheta,Double_t &phi,Float_t sat_ra[],Float_t sat_dec[],Double_t vector_in[]) {
  Float_t ux1[3];
  Float_t uy1[3];
  Float_t uz1[3];
  Float_t rax = sat_ra[0];
  Float_t ray = sat_ra[1];
  Float_t raz = sat_ra[2];
  Float_t decx = sat_dec[0];
  Float_t decy = sat_dec[1];
  Float_t decz = sat_dec[2];

  Double_t tmp_local[3],sinphi,cosphi;
  
  ux1[0] = cos(decx)*cos(rax);
  ux1[1] = cos(decx)*sin(rax);
  ux1[2] = sin(decx);

  uy1[0] = cos(decy)*cos(ray);
  uy1[1] = cos(decy)*sin(ray);
  uy1[2] = sin(decy);

  uz1[0] = cos(decz)*cos(raz);
  uz1[1] = cos(decz)*sin(raz);
  uz1[2] = sin(decz);
  
  tmp_local[2]=((ux1[0]*vector_in[1]-ux1[1]*vector_in[0])*(ux1[2]*uy1[0]-ux1[0]*uy1[2])+(ux1[0]*vector_in[2]-ux1[2]*vector_in[0])*(ux1[0]*uy1[1]-ux1[1]*uy1[0]))/((ux1[0]*uz1[2]-ux1[2]*uz1[0])*(ux1[0]*uy1[1]-ux1[1]*uy1[0])-(ux1[1]*uz1[0]-ux1[0]*uz1[1])*(ux1[2]*uy1[0]-ux1[0]*uy1[2]));

  tmp_local[1]=(tmp_local[2]*(ux1[1]*uz1[0]-ux1[0]*uz1[1])+ux1[0]*vector_in[1]-ux1[1]*vector_in[0])/((ux1[0]*uy1[1])-(ux1[1]*uy1[0]));

  tmp_local[0]=(1/ux1[0])*(vector_in[0]-tmp_local[1]*uy1[0]-tmp_local[2]*uz1[0]);

  costheta=tmp_local[2];
  sinphi=tmp_local[1]/sin(acos(costheta));
  cosphi=tmp_local[0]/sin(acos(costheta));

  if(sinphi>=0 && cosphi>=0)
    phi=acos(cosphi);
  else if(sinphi>0 && cosphi<0)
    phi=acos(cosphi);
  else if(sinphi<0 && cosphi<0) {
    phi=acos(cosphi);
    phi+=2*(TMath::Pi()-phi);
  }
  else {
    phi=acos(cosphi);
    phi=2*TMath::Pi()-phi;
  }
}
