#include <bitset>

#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedTowerCorrelator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

void l1tp2::DigitizedTowerCorrelator::print(const std::string location) {
  std::cout << "[Info:] DigitizedTowerCorrelator " << location << ": " << std::bitset<16>{data()} << std::endl;
}

