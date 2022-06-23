#include "Phase2L1RCT.h"
#include "Phase2L1GCT.h"

#ifndef _PHASE_2_L1_GCT_ALGO_H_
#define _PHASE_2_L1_GCT_ALGO_H_

// Get GCT cluster c's iEta (global iEta convention).
// Use with getEta_fromCrystaliEta from Phase2L1RCT.h to convert from GCT cluster
// to real eta.
int getCluster_global_iEta(unsigned int nGCTCard, GCTcluster_t c) {

  // First get the "iEta/iPhi" in the GCT card. i.e. in the diagram where the barrel
  // is split up into three GCT cards, (iEta, iPhi) = (0, 0) is the top left corner
  // of the GCT card. 
  int iEta_in_gctCard;

  // towEtaNeg = 1 (true) for negative eta
  if (c.towEtaNeg) {
    // c.towEta and c.crEta count outwards from the real eta = 0 center line, so for negative
    // eta, to convert to the barrel diagram global iEta (global iEta = 0 from LHS of page),
    // do (17*5 - 1) minus the GCT value.
    // e.g. If in GCT, a negative card's cluster had iEta = 84, this would be global iEta = 0.
    iEta_in_gctCard = ((17 * 5 - 1) - ((c.towEta * 5) + c.crEta));
  }
  else {
    // c.towEta and c.crEta count outwards from the real eta = 0 center line, so for positive
    // eta we need to add the 17*5 offset so that positive eta 0+epsilon starts at 17*5.
    // e.g. If in GCT, a positive card's cluster had iEta = 0, this would be global iEta = 85.
    // e.g. If in GCT, a positive card's cluster had iEta = 84, this would be global iEta = 169.
    iEta_in_gctCard = ((17 * 5) + ((c.towEta * 5) + c.crEta));
  }

  // Last, convert to the global iEta/iPhi in the barrel region. For eta there is nothing to
  // change.
  int iEta_in_barrel = iEta_in_gctCard;
  
  return iEta_in_barrel;
}

// Get GCT cluster c's iPhi (global convention).
// Use with getPhi_fromCrystaliPhi from Phase2L1RCT.h to convert from GCT cluster to real phi.
int getCluster_global_iPhi(unsigned int nGCTCard, GCTcluster_t c) {
  
  // First get the "iEta/iPhi" in the GCT card. i.e. in the diagram where the barrel
  // is split up into three GCT cards, (iEta, iPhi) = (0, 0) is the top left corner
  // of the GCT card.  

  // Luckily, in the GCT algo convention and the global convention, iPhi always increases from the
  // top of the page to the bottom of page in the barrel diagram.
  int iPhi_in_gctCard = ((c.towPhi * 5) + c.crPhi);

  // Last, convert to the global iEta/iPhi in the barrel region. For phi, we need to add the offset
  // of the GCT card in the phi direction, and modulo with the total number of crystals in the barrel
  // in the phi direction, since it wraps around.
  assert(nGCTCard <= 2);  // Make sure the card number is valid
  int iPhi_card_offset;
  if      (nGCTCard == 0) iPhi_card_offset = 20 * 5;  // tower #20, and five crystals per tower
  else if (nGCTCard == 1) iPhi_card_offset = 44 * 5;   
  else if (nGCTCard == 2) iPhi_card_offset = 68 * 5;

  int iPhi_in_barrel = (iPhi_card_offset + iPhi_in_gctCard) % (72 * 5); // detector wraps around in phi

  return iPhi_in_barrel;
}


GCTcard_t getClustersCombined(const GCTcard_t& GCTcard){

  GCTcard_t GCTcombinedClusters ;

  // first we initialize the output == input 
  //
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	
	// et2x5, et5x5, is_ss, and is_looseTkss remain the same whether or not the cluster is stitched across GCT card boundaries
	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ; 

	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;

	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;

	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss ;
      }}}

  // we will store new et in the GCTcombinedClusters, 0'ing lower clusters after stiching, dont need to care about other variables they stay the 
  // same as input for now at least
  // we combine even phi boudaries positive eta, when combined the lowest et is set to 0
  
  for(int i=0; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	ap_uint<15> phi1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	// if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4){
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
	      if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 ) {
		ap_uint<15> dPhi ; dPhi=((5 - phi1) + phi2) ;    
	        ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		// if( dEta < 2 ) {
		if( (dPhi <= 5) && (dEta < 2) ) {  
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
	// if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4) {
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3){ 
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if(GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
	      if(GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0) {  
		ap_uint<15> dPhi ; dPhi=((5 - phi1) + phi2) ; 
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( (dPhi <= 5) && (dEta < 2) ) {         
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
	// if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
	if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
	      if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3) {
		ap_uint<15> dPhi ; dPhi=((5 - phi2) + phi1) ;  // reversed for negative eta
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( (dPhi <= 5) && (dEta < 2) ) {    
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
	    
	// if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0){ 
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if( GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
	      if( GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 ) {
		ap_uint<15> dPhi ; dPhi=((5 - phi2) + phi1) ; // reversed compared to positive eta
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( (dPhi <= 5) && (dEta < 2) ) {        
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
	
	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  ;

	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  ;

	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  ;
	
	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  ;
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
	// positive eta: cluster eta/phi info
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et    = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et2x5 = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et5x5 = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].is_ss = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].is_looseTkss = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss;

	// positive eta: cluster eta/phi info
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towEtaNeg  = 0 ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towEta  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towPhi  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi+i*4  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].crEta  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta  ;
	GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].crPhi  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi  ;

	
	// negative eta: cluster energy info
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et    = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et2x5 = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et5x5 = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].is_ss = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;
	GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].is_looseTkss = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss;
	
	// negative eta: cluster eta/phi info
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

// algo_top: First two arguments are the same as in the original firmware.
// nGCTCard is 0, 1, or 2 (needed for getting the cluster real eta/phis for CMSSW collections).
// gctClusters is the CMSSW-style collection of clusters, to be used in the analyzer.
// gctTowers is the CMSSW-style collection of towers, to be used in the analyzer.

void algo_top(const GCTcard_t& GCTcard, GCTtoCorr_t& GCTtoCorr,
	      unsigned int nGCTCard,
	      std::unique_ptr<l1tp2::CaloCrystalClusterCollection> const& gctClusters,
	      std::unique_ptr<l1tp2::CaloTowerCollection> const& gctTowers) {
  
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
  // For the full GCT card (eight RCT cards wide in phi), skip the overlap region, i.e. SKIP i = 0, 1, 2, 3,
  // and i = 28, 29, 30, 31.
  for(int i=4; i<(N_GCTPOSITIVE_FIBERS-N_RCTGCT_FIBERS); i++){
    for(int k=0; k<N_GCTCLUSTERS_FIBER; k++){
      /* std::cout << "Accessing positive eta: GCTCorrfiber " << i-4  */
      /* 		<< " , GCTclusters " << k  */
      /* 		<< " , energy " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et << std::endl; */
      // Comments from Stephanie:
      // Use GCTinternal for GCT clusters. The indexing for towers and crystal eta/phi is that
      // etas of towers go from 0-16 where iEta = 0 is real eta = 0 (indexes increase towards
      // larger abs(eta). Phis of towers go from 0-3 where iPhi = 0 is the leftmost (if you look
      // at a 'sideways' diagram of the GCT card like on the TWiki). 
      // 
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et   ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towEtaNeg  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg  ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi-4 ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].crEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].crPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et2x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5 ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et5x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5 ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].is_ss  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].is_looseTkss = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss ; // new
      // Get the real eta, phi using two helper functions
      int globaliEta = getCluster_global_iEta(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      int globaliPhi = getCluster_global_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      float realEta = getEta_fromCrystaliEta(globaliEta);
      float realPhi = getPhi_fromCrystaliPhi(globaliPhi);
      
      reco::Candidate::PolarLorentzVector p4cluster(GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0,
						    realEta,
						    realPhi,
						    0.);
      l1tp2::CaloCrystalCluster cluster(p4cluster, 
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, // conver to float
					0,  // float h over e                              
                                        0,  // float iso                                       
                                        0,  // DetId seedCrystal                              
                                        0,  // puCorrPt                                           
                                        0,  // 0, 1, or 2 (as computed in firmware)                
                                        0,  // et2x2 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5/8.0,  // et2x5 (as computed in firmware, save float)           
                                        0,  // et3x5 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5/8.0,   // et5x5 (as computed in firmware, save float)  
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss,  // standalone WP
					false, // electronWP98: not computed
					false, // photonWP80: not computed
					false, // electronWP90: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss, // looseL1TkMatchWP
					false  // stage2effMatch: not computed
                                        );
      if (cluster.pt() > 0.0) {
	gctClusters->push_back(cluster);
	std::cout << "--- cluster pT, global iEta, iPhi and real eta, phi: "
		  << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0 
		  << ", ("
		  << globaliEta 
		  << ", "
	          << globaliPhi
		  << "), ("
		  << realEta 
		  << "," 
		  << realPhi
		  << ")"
	          << std::endl;
	std::cout << "    with the GCTinternal values: " 
		  << "towEtaNeg: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg << ", "
		  << "towEta: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta    << ", "
		  << "towPhi: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi    << ", "
		  << "crEta: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta   << ", "
		  << "crPhi: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi   << ", "
		  << "et2x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5/8.0 << ", "
		  << "et5x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5/8.0 << ", "
		  << "is_ss: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss << ", "
		  << "is_looseTkss" << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss << std::endl;
	
	
      }

    }
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      /* std::cout<< "Accessing positive eta: GCTCorrfiber " << i-4 */
      /* 	       << " , GCTtowers " << k  */
      /* 	       << " , energy " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].et << std::endl; */
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k].et  = GCTinternal.GCTCorrfiber[i].GCTtowers[k].et ;
    }
  }
  // In negative eta, skip the overlap region, i.e. SKIP i = 32, 33, 34, 35, and 61, 62, 63, 64.
  for(int i=(N_GCTPOSITIVE_FIBERS+N_RCTGCT_FIBERS); i<(N_GCTINTERNAL_FIBERS-N_RCTGCT_FIBERS); i++){
    for(int k=0; k<N_GCTCLUSTERS_FIBER; k++){
      /* std::cout << "Accessing negative eta: GCTCorrfiber " << i-12 */
      /*           << " , GCTclusters " << k  */
      /*           << " , energy " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et << std::endl; */
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].et  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et   ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towEtaNeg  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg  ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi-4 ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].crEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].crPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].et2x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5 ; // new                                                           
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].et5x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5 ; // new                                                             
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].is_ss  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss ; // new                                                 
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].is_looseTkss = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss ; // new    

      // Get the real eta, phi using two helper functions
      int globaliEta = getCluster_global_iEta(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      int globaliPhi = getCluster_global_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      float realEta = getEta_fromCrystaliEta(globaliEta);
      float realPhi = getPhi_fromCrystaliPhi(globaliPhi);
      
      reco::Candidate::PolarLorentzVector p4cluster(GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0,
						    realEta,
						    realPhi,
						    0.);
      l1tp2::CaloCrystalCluster cluster(p4cluster, 
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, // conver to float
					0,  // float h over e                              
                                        0,  // float iso                                       
                                        0,  // DetId seedCrystal                              
                                        0,  // puCorrPt                                           
                                        0,  // 0, 1, or 2 (as computed in firmware)                
                                        0,  // et2x2 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5/8.0,  // et2x5 (as computed in firmware, save float)           
                                        0,  // et3x5 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5/8.0,   // et5x5 (as computed in firmware, save float)  
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss,  // standalone WP
					false, // electronWP98: not computed
					false, // photonWP80: not computed
					false, // electronWP90: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss, // looseL1TkMatchWP
					false  // stage2effMatch: not computed
                                        );
      if (cluster.pt() > 0.0) {
	gctClusters->push_back(cluster);
	std::cout << "--- cluster pT, global iEta, iPhi and real eta, phi: "
		  << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0 
		  << ", ("
		  << globaliEta 
		  << ", "
		  << globaliPhi
		  << "), ("
		  << realEta 
		  << "," 
		  << realPhi
		  << ")"
		  << std::endl;
	std::cout << "    with the GCTinternal values: "
                  << "towEtaNeg: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg << ", "
                  << "towEta: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta    << ", "
                  << "towPhi: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi    << ", "
                  << "crEta: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta   << ", "
                  << "crPhi: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi   << ", "
		  << "et2x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5 << ", "
                  << "et5x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5 << ", "
                  << "is_ss: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss << ", "
                  << "is_looseTkss: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss << std::endl;
	
      }
      
    }
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      /* std::cout<< "Accessing positive eta: GCTCorrfiber " << i-4 */
      /* 	       << " , GCTtowers " << k */
      /* 	       << " , energy " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].et << std::endl; */
      GCTtoCorr.GCTCorrfiber[i-12].GCTtowers[k].et  = GCTinternal.GCTCorrfiber[i].GCTtowers[k].et ;
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

#endif
