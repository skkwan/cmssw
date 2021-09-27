#include "Phase2L1GCT.h"


GCTcard_t getClustersCombined(const GCTcard_t& GCTcard){

  GCTcard_t GCTcombinedClusters ;

  // first we initialize the output == input 
  //
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
      }}}

  // we will store new et in the GCTcombinedClusters, 0'ing lower clusters after stiching, dont need to care about other variables they stay the 
  // same as input for now at least
  // we combine even phi boudaries positive eta, when combined the lowest et is set to 0
  
  for(int i=0; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	ap_uint<15> phi1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4){
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( dEta < 2 ) {
		  ap_uint<12> one = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
		  ap_uint<12> two = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
		  if (one > two){ 
		    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
		    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
		  }
		  else {
		    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
		    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 

		  }
		}}
	    }}
	}
      }
    }
  }

  // now we combine odd phi boundaries positive eta
    
  for(int i=1; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	ap_uint<15> phi1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      if(GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( dEta <2 ) {
		  ap_uint<12> one = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
		  ap_uint<12> two = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
		  if (one > two){ 
		    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
		    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
		  }
		  else {
		    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
		    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
		  }
		}}
	    }}
	}
      }
    }
  }
    
  // repeat above steps for negative eta
  // even phi boundaries
    
  for(int i=0; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> eta1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	ap_uint<15> phi1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( dEta < 2 ) {
		  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
		  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
		  if (one > two){ 
		    GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
		    GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
		  }
		  else {
		    GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
		    GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
		  }
		}}
	    }}
	}
      }
    }
  }
    
  // now we combine odd boundaries 
    
  for(int i=1; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> eta1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	ap_uint<15> phi1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	    
	if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      if( GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if(dEta < 2 ) {
		  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
		  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
		  if (one > two){ 
		    GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
		    GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
		  }
		  else {
		    GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
		    GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
		  }
		}}
	    }}
	}
      }
    }
  }
  //
    
  // we need to store what we did before we start phi stiching
  //
  GCTcard_t GCTout ;
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
      }}}
          
  // now we combine eta boundaries, just positive and negative eta 
     
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> phi1 = (i*4+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi)*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta == 0 && eta1 == 0 ) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> phi2 = (i*4+(3-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi))*5+(4-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi) ;
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      if( GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].towEta == 0 && eta2 == 0 ) {
		ap_uint<15> dPhi ; dPhi=(phi1 > phi2)?(phi1-phi2):(phi2-phi1) ;
		if( dPhi < 2 ) {
		  ap_uint<12> one = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
		  ap_uint<12> two = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
		  if (one > two){ 
		    GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
		    GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
		  }
		  else {
		    GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
		    GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
		  }
		}}
	    }}
	}
      }
    }
  }
  return GCTout ;
}


GCTinternal_t getClustersTowers(const GCTcard_t& GCTcard){

  GCTcard_t GCTcombinedClusters ;
  GCTinternal_t GCTout ;

  // here we will stitch the clusters in phi and eta
  GCTcombinedClusters = getClustersCombined(GCTcard) ;

  // create internal structure of GCT card
  // we start from RCT card 0 - it is overlap with other GCT card and fill structure that we will use to send data to Correlator
  // we only need to care about clusters et in combinrdClusters, since the rest remains unchanged wrt input, the cluster that we set to 0
  // remain in the data at the same place , it will just get 0 et now
  // we need to code Positive and Negative Eta differently !  For negative Eta link 0 for each RCT
  // region becomes 3 in GCT output, the RCT card is rotated around 0:0 point of the card 
  // First 16 fibers - positive Eta , second 16 - negative. Eta coded 0...16 and towEtaNeg = 0 or 1 for clusters ; 
  // Phi is coded 0...15 , in case if whole card 0...33 and subdevision 1/5 in crPhi and crEta 0...4 for 
  // position in tower
  //
  // towers are put in link starting from eta=0, the link number defines Eta negative or positive and Phi position of tower.
  //
  //
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towEtaNeg  = 0 ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towEta  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towPhi  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi+i*4  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].crEta  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].crPhi  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi  ;
        
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].towEtaNeg  = 1 ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].towEta  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towEta  ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].towPhi  = (3-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi)+i*4  ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].crEta  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crEta  ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].crPhi  = (4-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi)  ;
      }
      for(int k=0; k<N_RCTTOWERS_FIBER; k++){
	GCTout.GCTCorrfiber[i*4+j].GCTtowers[k].et  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTtowers[k].et  ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTtowers[k].et  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTtowers[k].et  ;
      }
    }}
  return GCTout ;
}

GCTintTowers_t  getFullTowers(const GCTinternal_t& GCTinternal) {
  GCTintTowers_t GCTintTowers;
  
  //-- positive eta
  //
  for(int i=0; i<N_GCTPOSITIVE_FIBERS; i=i+4){
    for(int i1=0; i1<4; i1++){
      for(int k=0; k<N_GCTTOWERS_FIBER; k++){
	ap_uint<15> phi = i+i1 ; 
	ap_uint<15> eta = N_GCTETA/2 + k ; 
	GCTintTowers.GCTtower[eta][phi].et  = GCTinternal.GCTCorrfiber[phi].GCTtowers[k].et ;
	for(int ic1=0; ic1<4; ic1++){
	  for(int jc=0; jc<N_GCTCLUSTERS_FIBER; jc++){
	    ap_uint<15> eta1 = N_GCTETA/2 + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towEta ; 
	    ap_uint<15> phi1 = GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towPhi ; 
	    if( eta == eta1 && phi == phi1) GCTintTowers.GCTtower[eta][phi].et  = GCTintTowers.GCTtower[eta][phi].et + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].et ;
	  }
	}
      }
    }
  }

  //-- negative eta
  //
  for(int i=N_GCTPOSITIVE_FIBERS; i<N_GCTINTERNAL_FIBERS; i=i+4){
    for(int i1=0; i1<4; i1++){
      for(int k=0; k<N_GCTTOWERS_FIBER; k++){
	ap_uint<15> eta = N_GCTETA/2 - k - 1 ; 
	ap_uint<15> phi = i+i1-N_GCTPOSITIVE_FIBERS ; 
	GCTintTowers.GCTtower[eta][phi].et  = GCTinternal.GCTCorrfiber[i+i1].GCTtowers[k].et ;
	for(int ic1=0; ic1<4; ic1++){
	  for(int jc=0; jc<N_GCTCLUSTERS_FIBER; jc++){
	    ap_uint<15> eta1 = N_GCTETA/2 - 1 - GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towEta ;
	    ap_uint<15> phi1 = GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towPhi ;
	    if( eta == eta1 && phi == phi1) GCTintTowers.GCTtower[eta][phi].et  = GCTintTowers.GCTtower[eta][phi].et + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].et ;
	  }
	}
      }
    }
  }
  
  return GCTintTowers ;
  //end
}

//void algo_top(const GCTcard_t& GCTcard, ap_uint<15> cluster[10]){
void algo_top(const GCTcard_t& GCTcard, GCTtoCorr_t& GCTtoCorr){

  GCTinternal_t GCTinternal ;

  // first we fill the GCT area with proper et and eta/phi; the eta is 0....16 and towEtaNeg = 0 or 1
  // Phi 0....15 since this is half GCT card and 0...3 is overlap currently on one side, crEta/Phi 0...4
  
  GCTinternal = getClustersTowers(GCTcard) ;
  
  
  GCTintTowers_t GCTintTowers;
  // here we combine towers and clusters to get full towers, and move to internal GCT card coordinate 0...33 in eta and 0...15 in phi,
  // should change to 0...31 for the whole card
  
  GCTintTowers = getFullTowers(GCTinternal) ;

  //ok now we have towers combined with clusters, now need to recalculate positions of all clusters in GCT coordinates
  
  // OUTPUT TO CORRELATOR !!!!
  // remove overlap region, we send out 24 fibers: 12 pos/12 neg
  // all phi now decreases by 4
  // this we send to Correlator
  for(int i=4; i<N_GCTPOSITIVE_FIBERS; i++){
    for(int k=0; k<N_GCTCLUSTERS_FIBER; k++){
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et   ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towEtaNeg  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg  ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi-4 ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].crEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].crPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi ;
    }
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k].et  = GCTinternal.GCTCorrfiber[i].GCTtowers[k].et ;
    }
  }
  for(int i=N_GCTPOSITIVE_FIBERS+4; i<N_GCTINTERNAL_FIBERS; i++){
    for(int k=0; k<N_GCTCLUSTERS_FIBER; k++){
      GCTtoCorr.GCTCorrfiber[i-8].GCTclusters[k].et  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et   ;
      GCTtoCorr.GCTCorrfiber[i-8].GCTclusters[k].towEtaNeg  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg  ;
      GCTtoCorr.GCTCorrfiber[i-8].GCTclusters[k].towEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta ;
      GCTtoCorr.GCTCorrfiber[i-8].GCTclusters[k].towPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi-4 ;
      GCTtoCorr.GCTCorrfiber[i-8].GCTclusters[k].crEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta ;
      GCTtoCorr.GCTCorrfiber[i-8].GCTclusters[k].crPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi ;
    }
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      GCTtoCorr.GCTCorrfiber[i-8].GCTtowers[k].et  = GCTinternal.GCTCorrfiber[i].GCTtowers[k].et ;
    }
  }
  
  // this for test of the pattern that is in the tb file
  /*
    
    cluster[0] = GCTtoCorr.GCTCorrfiber[0].GCTclusters[0].et  ;
    cluster[1] = GCTtoCorr.GCTCorrfiber[0].GCTclusters[0].towEta  ;
    cluster[2] = GCTtoCorr.GCTCorrfiber[0].GCTclusters[0].towPhi  ;
    cluster[3] = 0  ;
    //           cluster[1] = GCTtoCorr.GCTCorrfiber[4].GCTclusters[0].et  ;
    //           cluster[2] = GCTtoCorr.GCTCorrfiber[15].GCTclusters[0].et  ;
    //           cluster[3] = GCTtoCorr.GCTCorrfiber[19].GCTclusters[0].et  ;
    //           cluster[4] = GCTtoCorr.GCTCorrfiber[18].GCTclusters[1].et  ;
    
    cluster[4] = 0  ;
    cluster[5] = 0  ;

    cluster[6] = GCTtoCorr.GCTCorrfiber[3].GCTtowers[7].et  ;
    cluster[7] = GCTtoCorr.GCTCorrfiber[10].GCTtowers[5].et  ;
    cluster[8] = GCTtoCorr.GCTCorrfiber[16].GCTtowers[7].et  ;
    cluster[9] = GCTtoCorr.GCTCorrfiber[23].GCTtowers[5].et  ;
    
    cluster[3] = GCTintTowers.GCTtower[12][8].et  ;
    cluster[4] = GCTintTowers.GCTtower[24][7].et  ;
    cluster[5] = GCTintTowers.GCTtower[18][4].et  ;
  */
  //end
}
