#ifndef L1Trigger_L1TTrackMatch_L1TkHTEmulatorProducer_HH
#define L1Trigger_L1TTrackMatch_L1TkHTEmulatorProducer_HH

#include <ap_int.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1Trigger/interface/TkJetWord.h"

// Namespace that defines constants and types used by the HT Emulation

namespace l1thtemu {

  const unsigned int kValidSize{1};
  const unsigned int kPtSize{18};
  const unsigned int kPtIntSize{13};
  const unsigned int kUnassignedSize{64 - (kPtSize + kValidSize)};

  enum BitLocations {
    // The location of the least significant bit (LSB) and most significant bit (MSB) in the sum word for different fields
    kValidLSB = 0,
    kValidMSB = kValidLSB + kValidSize - 1,
    kPtLSB = kValidMSB + 1,
    kPtMSB = kPtLSB + kPtSize - 1,
    kUnassignedLSB = kPtMSB + 1,
    kUnassignedMSB = kUnassignedLSB + kUnassignedSize - 1,
  };

  const float kMaxHT{8192};  // 8.192 TeV

  typedef ap_ufixed<kPtSize, kPtIntSize> ht_t;

  const unsigned int kHTBins = 1 << kPtSize;

  const double kStepPt{0.25};

  const double kStepHT = (l1thtemu::kMaxHT / l1thtemu::kHTBins);

  template <typename T>
  T digitizeSignedValue(double value, unsigned int nBits, double lsb) {
    T digitized_value = std::floor(std::abs(value) / lsb);
    T digitized_maximum = (1 << (nBits - 1)) - 1;  // The remove 1 bit from nBits to account for the sign
    if (digitized_value > digitized_maximum)
      digitized_value = digitized_maximum;
    if (value < 0)
      digitized_value = (1 << nBits) - digitized_value;  // two's complement encoding
    return digitized_value;
  }

  struct Ht {
    ht_t Et;
  };

}  // namespace l1thtemu
#endif
