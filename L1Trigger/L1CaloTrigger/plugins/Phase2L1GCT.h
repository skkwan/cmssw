#ifndef _PHASE_2_L1_GCT_H_
#define _PHASE_2_L1_GCT_H_

#include <iostream>
#include <ap_int.h>

static constexpr int N_RCTCARDS_PHI = 8;
static constexpr int N_RCTGCT_FIBERS = 4; 
static constexpr int N_RCTTOWERS_FIBER = 17;
static constexpr int N_RCTCLUSTERS_FIBER = 2;

static constexpr int N_GCTCARDS = 3;
static constexpr int N_GCTCORR_FIBERS = 48;
static constexpr int N_GCTTOWERS_FIBER = 17; 
static constexpr int N_GCTCLUSTERS_FIBER = 2; 

static constexpr int N_GCTINTERNAL_FIBERS = 32;
static constexpr int N_GCTPOSITIVE_FIBERS = 32;
static constexpr int N_GCTETA = 34;
static constexpr int N_GCTPHI = 32;

typedef ap_uint<5> loop;
using namespace std;

typedef struct {
  ap_uint<12> et ;
  ap_uint<5> towEta ;
  ap_uint<2> towPhi ;
  ap_uint<3> crEta ;
  ap_uint<3> crPhi ;
} RCTcluster_t ;

typedef struct {
  ap_uint<12> et ;
} RCTtower_t ;

typedef struct {
  RCTtower_t RCTtowers[N_RCTTOWERS_FIBER] ;
  RCTcluster_t RCTclusters[N_RCTCLUSTERS_FIBER] ;
} RCTtoGCTfiber_t ;

typedef struct {
  RCTtoGCTfiber_t RCTtoGCTfiber[N_RCTGCT_FIBERS] ;
} RCTcard_t ;

typedef struct {
  RCTcard_t RCTcardEtaPos[N_RCTCARDS_PHI] ;
  RCTcard_t RCTcardEtaNeg[N_RCTCARDS_PHI] ;
} GCTcard_t ;

typedef struct {
  ap_uint<12> et ;
  ap_uint<1> towEtaNeg ;
  ap_uint<6> towEta ;
  ap_uint<7> towPhi ;
  ap_uint<3> crEta ;
  ap_uint<3> crPhi ;
} GCTcluster_t ;

typedef struct {
  ap_uint<12> et ;
} GCTtower_t ;

typedef struct {
  GCTtower_t GCTtowers[N_GCTTOWERS_FIBER] ;
  GCTcluster_t GCTclusters[N_GCTCLUSTERS_FIBER] ;
} GCTCorrfiber_t ;

typedef struct {
  GCTCorrfiber_t GCTCorrfiber[N_GCTCORR_FIBERS] ;
} GCTtoCorr_t ;

typedef struct {
  GCTCorrfiber_t GCTCorrfiber[N_GCTINTERNAL_FIBERS] ;
} GCTinternal_t ;


typedef struct {
  GCTtower_t GCTtower[N_GCTETA][N_GCTPHI] ;
} GCTintTowers_t ;

/* For each GCT card (3 of them in total, for barrel + endcap), list the sixteen                
 * RCT cards that fall in them. The first eight are in positive eta, the next                   
 * eight are in negative eta (see figure of one GCT card). The RCT cards themselves             
 * run from 0 to 35 (see RCT figure).                                                          
 * Hard-coded because the GCT cards wrap around the barrel region.                            
 * Used only to convert the RCT emulator outputs to the GCT emulator inputs.                   
 */
static const unsigned int GCTcardtoRCTcardnumber[N_GCTCARDS][N_RCTCARDS_PHI * 2] =
  { // GCT Card 0                                                                            
    {11, 13, 15, 17, 19, 21, 23, 25,
     10, 12, 14, 16, 18, 20, 22, 24},

    // GCT Card 1                                                                            
    {23, 25, 27, 29, 31, 33, 35,  1,
     22, 24, 26, 28, 30, 32, 34,  0},

    // GCT Card 2                                                                            
    {35,  1,  3,  5,  7,  9, 11, 13,
     34,  0,  2,  4,  6,  8, 10, 12}
  };

/*
 * Initialize a GCT card so that its clusters and towers have zero energy.
 */
void initializeGCTCard(GCTcard_t gctCard){
  for (int i = 0; i < (N_RCTCARDS_PHI); i++) {  // loop from i = 0 to i = 8 and do pos and neg eta at the same time
    for (int iLink = 0; iLink < N_RCTGCT_FIBERS; ++iLink) {  // four links per card
      for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; ++iTower) {  // 17 towers per link
	RCTtower_t t0;
	t0.et = 0;
	gctCard.RCTcardEtaPos[i].RCTtoGCTfiber[iLink].RCTtowers[iTower] = t0;
	gctCard.RCTcardEtaNeg[i].RCTtoGCTfiber[iLink].RCTtowers[iTower] = t0;
      }
      
      for (int iCluster = 0; iCluster < N_GCTCLUSTERS_FIBER; ++iCluster) { // 2 clusters per link
	RCTcluster_t c0;
	c0.et     = 0;
	c0.towEta = 0;
	c0.towPhi = 0;
	c0.crEta  = 0;
	c0.crPhi  = 0;
	gctCard.RCTcardEtaPos[i].RCTtoGCTfiber[iLink].RCTclusters[iCluster] = c0;
	gctCard.RCTcardEtaNeg[i].RCTtoGCTfiber[iLink].RCTclusters[iCluster] = c0;
      }
    }
  }
}

//void algo_top(const GCTcard_t& GCTcard, ap_uint<15> cluster[10]);
void algo_top(const GCTcard_t& GCTcard, GCTtoCorr_t& GCTtoCorr) ;

GCTinternal_t getClustersTowers(const GCTcard_t& GCTcard) ;

GCTcard_t getClustersCombined(const GCTcard_t& GCTcard) ; 

GCTintTowers_t  getFullTowers(const GCTinternal_t& GCTinternal) ;

#endif
