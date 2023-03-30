#include <bitset>

#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterCorrelator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

void l1tp2::DigitizedClusterCorrelator::print(const std::string location) {
  //edm::LogInfo("DigitizedClusterCorrelator") << "Info: DigitizedClusterCorrelator " << location << ": " << (int) data() << std::endl;
  std::cout << "[Info:] DigitizedClusterCorrelator " << location << ": " << std::bitset<64>{data()} << std::endl;
}
