#include <bitset>

#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedTowerGT.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

void l1tp2::DigitizedTowerGT::print(const std::string location) {
  std::cout << "[Info:] DigitizedTowerGT " << location << ": " << std::bitset<16>{data()} << std::endl;
}

