#include <bitset>

#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterGT.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

void l1tp2::DigitizedClusterGT::print(const std::string location) {
  std::cout << "[Info:] DigitizedClusterGT " << location << ": " << std::bitset<64>{data()} << std::endl;
}
