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

  // This needs to mimic the firmware object

  const unsigned int kValidSize{1};
  const unsigned int kVectorSumSize{16}; // unused: vector sum is ap_ufixed(16, 11)
  const unsigned int kVectorSumIntSize{11}; // unused: vector sum is ap_ufixed(16, 11)
  const unsigned int kVectorSumPhiSize{13}; // unused: vector sum phi is 13 bits
  const unsigned int kScalarSumHTSize{18};
  const unsigned int kScalarSumHTIntSize{13};
  const unsigned int kUnassignedSize{64 - (kScalarSumHTSize + kVectorSumSize + kVectorSumPhiSize + kValidSize)};

  enum BitLocations {
    // The location of the least significant bit (LSB) and most significant bit (MSB) in the sum word for different fields
    kValidLSB = 0,
    kValidMSB = kValidLSB + kValidSize - 1,
    kVectorSumLSB = kValidMSB + 1,
    kVectorSumMSB = kVectorSumLSB + kVectorSumSize - 1,
    kVectorSumPhiLSB = kVectorSumMSB + 1,
    kVectorSumPhiMSB = kVectorSumPhiLSB + kVectorSumPhiSize - 1,
    kScalarSumHTLSB = kVectorSumPhiMSB + 1,
    kScalarSumHTMSB = kScalarSumHTLSB + kScalarSumHTSize - 1,
    kUnassignedLSB = kScalarSumHTMSB + 1,
    kUnassignedMSB = kUnassignedLSB + kUnassignedSize - 1,
  };

  const float kMaxHT{8192};  // 8.192 TeV

  typedef ap_ufixed<kScalarSumHTSize, kScalarSumHTIntSize> ht_t;

  const unsigned int kHTBins = 1 << kScalarSumHTSize;

  const double kStepPt = 0.03125; // jet pT
  const unsigned int kPtSize{16}; // jet pT int size ap_ufixed(16, 11)

  const double kStepHT = (l1thtemu::kMaxHT / l1thtemu::kHTBins); // (8192 / (1<<18)) = 0.03125

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
