#ifndef _PHASE_2_L1_GCT_H_
#define _PHASE_2_L1_GCT_H_

#include <iostream>
#include <ap_int.h>

#define N_RCTCARDS_PHI 8 
#define N_RCTGCT_FIBERS 4 
#define N_RCTTOWERS_FIBER 17 
#define N_RCTCLUSTERS_FIBER 2 

#define N_GCTCARDS 3
#define N_GCTCORR_FIBERS 48
#define N_GCTTOWERS_FIBER 17 
#define N_GCTCLUSTERS_FIBER 2 

#define N_GCTINTERNAL_FIBERS 32
#define N_GCTPOSITIVE_FIBERS 32
#define N_GCTETA 34
#define N_GCTPHI 32

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


//void algo_top(const GCTcard_t& GCTcard, ap_uint<15> cluster[10]);
void algo_top(const GCTcard_t& GCTcard, GCTtoCorr_t& GCTtoCorr) ;

GCTinternal_t getClustersTowers(const GCTcard_t& GCTcard) ;

GCTcard_t getClustersCombined(const GCTcard_t& GCTcard) ; 

GCTintTowers_t  getFullTowers(const GCTinternal_t& GCTinternal) ;

#endif
