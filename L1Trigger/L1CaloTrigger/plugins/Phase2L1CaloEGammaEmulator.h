//------------------------------------
// EGammaCrystalsProducer.h 
//------------------------------------

#include <fstream>
#include <stdio.h>
#include <iostream>



#ifndef EGAMMA_CRYSTALS_PRODUCER_H
#define EGAMMA_CRYSTALS_PRODUCER_H

static constexpr int dimPhi = 4;
static constexpr int dimEta = 17;
static constexpr int dimCard = 36;

//--------------------------------------------------------// 

// Helper function to print a 4x17x36 int array to a ostream f.                                                                

inline void printL1ArrayInt(ofstream& f, int array[dimPhi][dimEta][dimCard],
                     string desc = "") {

  for (int kk = 0; kk < dimCard; kk++) {
    f << "[CARD " << kk << "]: " << desc << std::endl;
    for (int ii = 0; ii < dimPhi; ii++) {
      for (int jj = 0; jj < dimEta; jj++) {
        f << array[ii][jj][kk] << "\t";
      }
      f << std::endl;
    }
  }

}

//--------------------------------------------------------//                                                                  

// Helper function to print a 4x17x36 float array to a ostream f.                                                                    
          
inline void printL1ArrayFloat(ofstream& f, float array[dimPhi][dimEta][dimCard],
                       string desc = "") {

  for (int kk = 0; kk < dimCard; kk++) {
    f << "[CARD " << kk << "]: " << desc << std::endl;
    for (int ii = 0; ii < dimPhi; ii++) {
      for (int jj = 0; jj < dimEta; jj++) {
        f << array[ii][jj][kk] << "\t";
      }
      f << std::endl;
    }
  }
}


#endif
