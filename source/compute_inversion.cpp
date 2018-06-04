
#include "MyHead.h"

///////////////////////////////////////////// Math custom functions

static Double_t acot(Double_t value) { return atan(1/value); }

//////////////////////////////////////////////////////////////////////////////

void sky_backtrack_rndm(Float_t sat_ra[],Float_t sat_dec[],TH2D *evDist_border,TH2D &inverse,Bool_t &first_call,Int_t &n_peacks,Int_t &n_valleys,std::ofstream &log_file,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b) {
    
    Double_t costheta=-999,phi=-999;
    
    Double_t max_b=0,min_b=0,max_l=0,min_l=0;
    Double_t Nmax_b=0,Nmin_b=0,Pmax_b=0,Pmin_b=0;
    Bool_t InsideFrame = true;
    Double_t choosen_l=0,choosen_b=0;
    
    TRandom3 r_gen(random_seed);
    
    //////////////////////////////////////////////
    
    for(Int_t idx_ev=0; idx_ev<sky_events; idx_ev++) {
        if(!POI_discrete_plotting(sat_ra,sat_dec,evDist_border,peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,first_call,n_peacks,n_valleys)) {
            std::cout<<"\nJumped event "<<idx_ev<<" ! Event Distribution pattern ricognition error\n";
            log_file <<"\nJumped event "<<idx_ev<<" ! Event Distribution pattern ricognition error\n";
            continue;
        }
        if(points_are_inside(Gpeack_l,Gvalley_l,n_peacks,n_valleys)==true) {
            obtain_galactic_max_min(Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,n_peacks,n_valleys,true,max_b,min_b,max_l,min_l,Nmax_b,Nmin_b,Pmax_b,Pmin_b);
            InsideFrame=true;
        }
        else {
            obtain_galactic_max_min(Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,n_peacks,n_valleys,false,max_b,min_b,max_l,min_l,Nmax_b,Nmin_b,Pmax_b,Pmin_b);
            InsideFrame=false;
        }
        if(InsideFrame) {
            choosen_l=r_gen.Uniform(min_l,max_l);
            choosen_b=r_gen.Uniform(min_b,max_b);
        }
        R_invert_map(costheta,phi,choosen_l,choosen_b,sat_ra,sat_dec);
        if(outside_evDist(costheta,phi,evDist_border)==true)
            idx_ev--;
        else
            inverse.Fill(costheta,phi);
    }
}

bool POI_discrete_plotting(Float_t sat_ra[],Float_t sat_dec[],TH2D* evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Bool_t &first_call,Int_t &n_peacks,Int_t &n_valleys) {
    
    Bool_t good_PR=true;   //Bool variable. True if the pattern recognition of the event distribution works, false in the other case.
    
    // This kind of study, at the countrary, should be done on all DAMPE acquisition period
    if(first_call) {
        //calculate_POI_stuff(sat_ra,sat_dec,evDist_border,peacks_costheta,peacks_phi,valley_costheta,valley_phi);
        get_relevant_evDist_points(evDist_border,peacks_costheta,valley_costheta,peacks_phi,valley_phi);
        n_peacks=peacks_costheta.size();
        n_valleys=valley_costheta.size();
        
        Gpeack_l.resize(n_peacks);
        Gpeack_b.resize(n_peacks);
        Gvalley_l.resize(n_valleys);
        Gvalley_b.resize(n_valleys);
        
        GetPOI(peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,sat_ra,sat_dec,good_PR,n_peacks,n_valleys);
        first_call=false;
    }
    
    GetPOI(peacks_costheta,peacks_phi,valley_costheta,valley_phi,Gpeack_l,Gpeack_b,Gvalley_l,Gvalley_b,sat_ra,sat_dec,good_PR,n_peacks,n_valleys);
    
    return good_PR;
    
}

void get_evDist_border(TH2D &evDist,TH2D &evDist_border) {
    
    for(Int_t y_bin=1; y_bin<=evDist.GetNbinsY(); y_bin++)
        for(Int_t x_bin=1; x_bin<=evDist.GetNbinsX(); x_bin++)
            if(evDist.GetBinContent(x_bin,y_bin)!=0) {
                evDist_border.SetBinContent(x_bin,y_bin,evDist.GetBinContent(x_bin,y_bin));
                break; //I found the first not-empty bin regarding a such phi (or y) value. Now I have to choose a new phi bin and search for the first not-empty costheta bin
            }
}

void get_relevant_evDist_points(TH2D* evDist_border,std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_phi) {
    
    Double_t l_binX,l_binY,costheta,phi;
    Double_t new_costheta,new_phi;
    Int_t bunch=get_bunch_dimension(evDist_border),n_bunch=evDist_border->GetNbinsY()/bunch;
    Bool_t first_b_bin,first_b=true;
    
    Double_t tmp_min_costheta1,tmp_max_costheta1,tmp_min_phi1,tmp_max_phi1;
    Double_t tmp_min_costheta2,tmp_max_costheta2,tmp_min_phi2,tmp_max_phi2;
    Double_t tmp_min_costheta3,tmp_max_costheta3,tmp_min_phi3,tmp_max_phi3;
    
    l_binX=(Double_t).5/evDist_border->GetNbinsX();
    l_binY=(Double_t)2*TMath::Pi()/evDist_border->GetNbinsY();
    
    for(Int_t idx_b=0; idx_b<n_bunch; idx_b++) {
        first_b_bin=true;
        //for(Int_t idx_bY=(((idx_b+3)*bunch)+1); idx_bY<=((idx_b+4)*bunch); idx_bY++) {
        for(Int_t idx_bY=((idx_b*bunch)+1); idx_bY<=((idx_b+1)*bunch); idx_bY++) {
            for(Int_t idx_bX=1; idx_bX<=evDist_border->GetNbinsX(); idx_bX++) {
                if(evDist_border->GetBinContent(idx_bX,idx_bY)!=0) {
                    costheta=0.5+.5*(idx_bX*l_binX+(idx_bX-1)*l_binX);
                    phi=.5*(idx_bY*l_binY+(idx_bY-1)*l_binY);
                    if(first_b_bin) {
                        tmp_min_costheta1=tmp_max_costheta1=costheta;
                        tmp_min_phi1=tmp_max_phi1=phi;
                        first_b_bin=false;
                    }
                    else {
                        if(costheta>tmp_max_costheta1) {
                            tmp_max_costheta1=costheta;
                            tmp_max_phi1=phi;
                        }
                        if(costheta<tmp_min_costheta1) {
                            tmp_min_costheta1=costheta;
                            tmp_min_phi1=phi;
                        }
                    }
                }
            }
        }
        first_b_bin=true;
        //for(Int_t idx_bY=(((idx_b+3)*bunch)+1); idx_bY<=((idx_b+4)*bunch); idx_bY++) {
        for(Int_t idx_bY=(((idx_b+1)*bunch)+1); idx_bY<=((idx_b+2)*bunch); idx_bY++) {
            for(Int_t idx_bX=1; idx_bX<=evDist_border->GetNbinsX(); idx_bX++) {
                if(evDist_border->GetBinContent(idx_bX,idx_bY)!=0) {
                    new_costheta=0.5+.5*(idx_bX*l_binX+(idx_bX-1)*l_binX);
                    new_phi=.5*(idx_bY*l_binY+(idx_bY-1)*l_binY);
                    if(first_b_bin) {
                        tmp_min_costheta3=tmp_max_costheta3=new_costheta;
                        tmp_min_phi3=tmp_max_phi3=new_phi;
                        first_b_bin=false;
                    }
                    else {
                        if(costheta>tmp_max_costheta3) {
                            tmp_max_costheta3=new_costheta;
                            tmp_max_phi3=new_phi;
                        }
                        if(costheta<tmp_min_costheta3) {
                            tmp_min_costheta3=new_costheta;
                            tmp_min_phi3=new_phi;
                        }
                    }
                }
            }
        }
        if(first_b)
            first_b=false;
        else {
            if(tmp_min_costheta1<tmp_min_costheta2 && tmp_min_costheta1<tmp_min_costheta3) {
                peacks_costheta.push_back(tmp_min_costheta1);
                peacks_phi.push_back(tmp_min_phi1);
            }
            if(tmp_max_costheta1>tmp_max_costheta2 && tmp_max_costheta1>tmp_max_costheta3) {
                valley_costheta.push_back(tmp_max_costheta1);
                valley_phi.push_back(tmp_max_phi1);
            }
        }
        tmp_min_costheta2=tmp_min_costheta1;
        tmp_max_costheta2=tmp_max_costheta1;
        tmp_min_phi2=tmp_min_phi1;
        tmp_max_phi2=tmp_max_phi1;
    }
    
}

Int_t get_bunch_dimension(TH2D* evDist_border) {
    Int_t bin_num=20,Ybins=evDist_border->GetNbinsY();
    Bool_t found=false;
    
    while(found==false)
        if((Ybins%bin_num)==0)
            found=true;
        else
            bin_num++;
    
    return bin_num;
}

void GetPOI(std::vector<Double_t> &peacks_costheta,std::vector<Double_t> &peacks_phi,std::vector<Double_t> &valley_costheta,std::vector<Double_t> &valley_phi,std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Float_t sat_ra[],Float_t sat_dec[],Bool_t &good_PR,Int_t n_peacks,Int_t n_valleys) {
    
    Double_t l=0,b=0;
    
    for(Int_t idx_p=0; idx_p<n_peacks; idx_p++) {
        R_from_local_to_galactic(peacks_costheta.at(idx_p),peacks_phi.at(idx_p),l,b,sat_ra,sat_dec);
        if(isnan(l)) {
            good_PR=false;
            break;
        }
        if(l>180.0)
            l-=360.0;
        Gpeack_b[idx_p]=b;
        Gpeack_l[idx_p]=l;
    }
    
    if(good_PR) {
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            R_from_local_to_galactic(valley_costheta.at(idx_v),valley_phi.at(idx_v),l,b,sat_ra,sat_dec);
            if(isnan(l)) {
                good_PR=false;
                break;
            }
            if(l>180.0)
                l-=360.0;
            Gvalley_b[idx_v]=b;
            Gvalley_l[idx_v]=l;
        }
    }
}

Double_t obtain_min_histo(TH1D *histo) {
    Double_t min_histo=-999,bin_lenght;
    
    bin_lenght=(Double_t).5/histo->GetNbinsX();
    
    for(Int_t idx_b=1; idx_b<=histo->GetNbinsX(); idx_b++)
        if(histo->GetBinContent(idx_b)!=0) {
            min_histo=0.5*(idx_b*bin_lenght+(idx_b-1)*bin_lenght);
            break; //the first not-empty bin has been found. No need to go ahead
        }
    
    return min_histo;
}


/////////////////////////////////////////////////////////////////////////////////////////



bool points_are_inside(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gvalley_l,Int_t n_peacks,Int_t n_valleys) {
    bool all_points_inside=true;
    
    /*
     for(Int_t idx_v=0; idx_v<5; idx_v++) {
     if(idx_v==0) {
     if(valley_glon[idx_v]>peacks_glon[0]) {
     all_points_inside=false;
     break;
     }
     }
     if(idx_v>0 && idx_v<3) {
     if(peacks_glon[idx_v]>peacks_glon[0] || valley_glon[idx_v]>peacks_glon[0]) {
     all_points_inside=false;
     break;
     }
     }
     else {
     if(peacks_glon[idx_v]>peacks_glon[0]) {
     all_points_inside=false;
     break;
     }
     }
     }
     */
    
    Double_t max=Gpeack_l[0],min=max;
    
    for(Int_t idx_v=1; idx_v<n_peacks; idx_v++) {
        if(Gpeack_l[idx_v]>max)
            max=Gpeack_l[idx_v];
        if(Gpeack_l[idx_v]<min)
            min=Gpeack_l[idx_v];
    }
    
    for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
        if(Gvalley_l[idx_v]>max)
            max=Gvalley_l[idx_v];
        if(Gvalley_l[idx_v]<min)
            min=Gvalley_l[idx_v];
    }
    
    if((max-min)>180.)
        all_points_inside=false;
    
    return all_points_inside;
}

void obtain_galactic_max_min(std::vector<Double_t> &Gpeack_l,std::vector<Double_t> &Gpeack_b,std::vector<Double_t> &Gvalley_l,std::vector<Double_t> &Gvalley_b,Int_t n_peacks,Int_t n_valleys,bool points_inside,Double_t &max_b,Double_t &min_b,Double_t &max_l,Double_t &min_l,Double_t &Nmax_b,Double_t &Nmin_b,Double_t &Pmax_b,Double_t &Pmin_b) {
    
    if(points_inside==true) {
        max_l=min_l=Gpeack_l[0];
        for(Int_t idx_v=1; idx_v<n_peacks; idx_v++) {
            if(Gpeack_l[idx_v]>max_l)
                max_l=Gpeack_l[idx_v];
            if(Gpeack_l[idx_v]<min_l)
                min_l=Gpeack_l[idx_v];
        }
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            if(Gvalley_l[idx_v]>max_l)
                max_l=Gvalley_l[idx_v];
            if(Gvalley_l[idx_v]<min_l)
                min_l=Gvalley_l[idx_v];
        }
        
        max_b=min_b=Gpeack_b[0];
        for(Int_t idx_v=1; idx_v<n_peacks; idx_v++) {
            if(Gpeack_b[idx_v]>max_b)
                max_b=Gpeack_b[idx_v];
            if(Gpeack_b[idx_v]<min_b)
                min_b=Gpeack_b[idx_v];
        }
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            if(Gvalley_b[idx_v]>max_b)
                max_b=Gpeack_b[idx_v];
            if(Gpeack_b[idx_v]<min_b)
                min_b=Gpeack_b[idx_v];
        }
    }
    else {
        max_l=-180;
        min_l=180;
        Nmax_b=Pmax_b=-90;
        Nmin_b=Pmin_b=90;
        
        for(Int_t idx_v=0; idx_v<n_peacks; idx_v++) {
            if(Gpeack_l[idx_v]<0) {
                if(Gpeack_l[idx_v]>max_l) {
                    max_l=Gpeack_l[idx_v];
                }
                if(Gpeack_b[idx_v]>Nmax_b) {
                    Nmax_b=Gpeack_b[idx_v];
                }
                if(Gpeack_b[idx_v]<Nmin_b) {
                    Nmin_b=Gpeack_b[idx_v];
                }
            }
            else {
                if(Gpeack_l[idx_v]<min_l) {
                    min_l=Gpeack_l[idx_v];
                }
                if(Gpeack_b[idx_v]>Pmax_b) {
                    Pmax_b=Gpeack_b[idx_v];
                }
                if(Gpeack_b[idx_v]<Pmin_b) {
                    Pmin_b=Gpeack_b[idx_v];
                }
            }
        }
        
        for(Int_t idx_v=0; idx_v<n_valleys; idx_v++) {
            if(Gvalley_l[idx_v]<0) {
                if(Gvalley_l[idx_v]>max_l) {
                    max_l=Gvalley_l[idx_v];
                }
                if(Gvalley_b[idx_v]>Nmax_b) {
                    Nmax_b=Gvalley_b[idx_v];
                }
                if(Gvalley_b[idx_v]<Nmin_b) {
                    Nmin_b=Gvalley_b[idx_v];
                }
            }
            else {
                if(Gvalley_l[idx_v]<min_l) {
                    min_l=Gvalley_l[idx_v];
                }
                if(Gvalley_b[idx_v]>Pmax_b) {
                    Pmax_b=Gvalley_b[idx_v];
                }
                if(Gvalley_b[idx_v]<Pmin_b) {
                    Pmin_b=Gvalley_b[idx_v];
                }
            }
        }
    }
}

void sky_backtrack(Float_t sat_ra[],Float_t sat_dec[],TH2D* evDist,TH2D *evDist_border,TH2D &direct,TH2D &inverse,TH1D &h_cos,TH1D &h_phi,TH1D &h_cos_LE,TH1D &h_phi_LE,Double_t &max_costheta,Double_t &min_costheta,Double_t &max_phi,Double_t &min_phi,Double_t tree_idx) {
  
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
    
        evDist->GetRandom2(costheta,phi);

        right_costheta=costheta;
        right_phi=phi;

        from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec,right_ra,right_dec,right_vector);  //Ok, now I have the coordinates into the galactic frame. I try to go back returning to the local frame.

        //////////// Control for "nan" longitude values values
        if(isnan(l)) {
            idx_t--;
            continue;
        }
      
        invert_map(costheta,phi,l,b,sat_ra,sat_dec,inv_ra,inv_dec,inv_vector);
        if(outside_evDist(costheta,phi,evDist_border)==true) {
            idx_t--;
            continue;
        }
          
        direct.Fill(right_costheta,right_phi);
        inverse.Fill(costheta,phi);
    
        diff_costheta=right_costheta-costheta;
        diff_phi=right_phi-phi;
      
        h_cos.Fill(diff_costheta);
        h_phi.Fill(diff_phi);
        h_cos_LE.Fill(diff_costheta);
        h_phi_LE.Fill(diff_phi);
      
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

    //obtain_costheta_phi(costheta,phi,sat_ra,sat_dec,vector_in);
    obtain_costheta_phi_ROOTf(costheta,phi,sat_ra,sat_dec,vector_in);
}

void R_invert_map(Double_t &costheta,Double_t &phi,Double_t l,Double_t b,Float_t sat_ra[],Float_t sat_dec[]) {
    Double_t ra=-999, dec=-999;
    AtPolarVect vector_out, vector_out_inv;
    AtVect vector_in;
    
    from_galactic_to_celestial(ra,dec,l,b);
    
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
    
    ///// Now obtain costheta and phi of the icnoming CR to the satellite
    
    //obtain_costheta_phi(costheta,phi,sat_ra,sat_dec,vector_in);
    obtain_costheta_phi_ROOTf(costheta,phi,sat_ra,sat_dec,vector_in);
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

void obtain_costheta_phi_ROOTf(Double_t &costheta,Double_t &phi,Float_t sat_ra[],Float_t sat_dec[],Double_t vector_in[]) {
    Float_t ux1[3];
    Float_t uy1[3];
    Float_t uz1[3];
    Float_t rax = sat_ra[0];
    Float_t ray = sat_ra[1];
    Float_t raz = sat_ra[2];
    Float_t decx = sat_dec[0];
    Float_t decy = sat_dec[1];
    Float_t decz = sat_dec[2];
    
    Double_t tmp_local[3],sinphi,cosphi,tmp_sum;
    TMatrixD Us(3,3),Us_invert(3,3);
    
    ux1[0] = cos(decx)*cos(rax);
    ux1[1] = cos(decx)*sin(rax);
    ux1[2] = sin(decx);
    
    uy1[0] = cos(decy)*cos(ray);
    uy1[1] = cos(decy)*sin(ray);
    uy1[2] = sin(decy);
    
    uz1[0] = cos(decz)*cos(raz);
    uz1[1] = cos(decz)*sin(raz);
    uz1[2] = sin(decz);
    
    //Fill matrix with Us already calculated
    
    Us(0,0)=ux1[0];
    Us(0,1)=ux1[1];
    Us(0,2)=ux1[2];
    
    Us(1,0)=uy1[0];
    Us(1,1)=uy1[1];
    Us(1,2)=uy1[2];
    
    Us(2,0)=uz1[0];
    Us(2,1)=uz1[1];
    Us(2,2)=uz1[2];
    
    ////// Invert Us
    
    Us_invert=Us.Invert();
    
    ////// Obtain tmp_local array
    
    for(Int_t idx_c=0; idx_c<3; idx_c++) {
        tmp_sum=0;
        for(Int_t idx_r=0; idx_r<3; idx_r++)
            tmp_sum+=vector_in[idx_r]*Us_invert(idx_r,idx_c);
        tmp_local[idx_c]=tmp_sum;
    }
    
    ///// Obtain costheta and phi in local reference frame !
    
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

void get_acceptance_border(TH2D *histo,TH2D* histo_border) {
    
    for(Int_t y_bin=1; y_bin<=histo->GetNbinsY(); y_bin++)
        for(Int_t x_bin=1; x_bin<=histo->GetNbinsX(); x_bin++)
            if(histo->GetBinContent(x_bin,y_bin)!=0) {
                histo_border->SetBinContent(x_bin,y_bin,histo->GetBinContent(x_bin,y_bin));
                break; //I found the first not-empty bin regarding a such phi (or y) value. Now I have to choose a new phi bin and search for the first not-empty costheta bin
            }
}

Bool_t outside_evDist(Double_t costheta,Double_t phi,TH2D* evDist_border) {
    Bool_t outside;
    Int_t Ybin;
    Double_t evDist_costheta=0,lXbin=0;
    
    TAxis *Yaxis = evDist_border->GetYaxis();
    TAxis *Xaxis = evDist_border->GetXaxis();
    Ybin=Yaxis->FindBin(phi);
    lXbin=.5/Xaxis->GetNbins();
    
    for(Int_t x_bin=1; x_bin<=evDist_border->GetNbinsX(); x_bin++)
        if(evDist_border->GetBinContent(x_bin,Ybin)!=0) {
            //evDist_costheta=x_bin*(.5/Xaxis->GetNbins());
            evDist_costheta = .5*((x_bin+1)*lXbin-x_bin*lXbin);
            break;
        }
    if(evDist_costheta>costheta)
        outside=true;
    else
        outside=false;
    
    return outside;
}
