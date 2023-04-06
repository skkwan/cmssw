#include <bitset>
#include <cassert>
#include <cmath>
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterGT.h"

/*
 * Helper for tests
 */
double round_to(double value, double precision = 1.0) { return std::round(value / precision) * precision; }

/*
 * Test basic constructor from digitized inputs
 */
int passItemizedDigitizedConstructorTest(void) {
  bool isFullyDigitizedInputs = true;
  ap_uint<1> isValid = 0x1;
  ap_uint<16> pt = 0x1;
  ap_uint<13> phi = 0x1;
  ap_uint<14> eta = 0x1;
  l1tp2::DigitizedClusterGT cGT = l1tp2::DigitizedClusterGT(isValid, pt, phi, eta, isFullyDigitizedInputs);
  return (cGT.passNullBitsCheck() && (cGT.pt() == 0x1) && (cGT.phi() == 0x1) && (cGT.eta() == 0x1));
}

/*
 * Positive saturated values constructor
 */
int passPositiveConstructorTest(void) {
  bool isValid = true;
  float pt_f = 2048;     // max value is 2047.97
  float phi_f = M_PI;    // phi = +180 should be 0b 1 1111 1111 1110 (LSB is the sign)
  float eta_f = 1.4841;  // becomes 1.4841 / (pi/ (2^12 - 1)) = 1934.49, truncates to 1934 = 0x78E = 0b11110001110
  // (contd.) finally bit shift left by 1, with the LSB = 0 since this is positive, which gives 0b111100011100 = 0xF1C

  l1tp2::DigitizedClusterGT cGT = l1tp2::DigitizedClusterGT(isValid, pt_f, phi_f, eta_f);
  return ((cGT.pt() == 0xFFFF) && (cGT.phi() == 0x1FFE) && (cGT.eta() == 0xF1C) && cGT.passNullBitsCheck());
}

/*
 * Negative saturated values constructor
 */
int passNegativeConstructorTest(void) {
  bool isValid = true;
  float pt_f = 2048;        // max value is 2047.97
  float phi_f = -1 * M_PI;  // phi = -180 should be 0b 1 1111 1111 1111 (LSB is the sign) = 0x1FFF
  float eta_f = -1.4841;    // becomes 1.4841 / (pi/ (2^12 - 1)) = 1934.49, truncates to 1934 = 0x78E = 0b11110001110
  // (contd.) finally bit shift left by 1, with the LSB = 1 since this is positive, which gives 0b111100011101 = 0xF1D

  l1tp2::DigitizedClusterGT cGT = l1tp2::DigitizedClusterGT(isValid, pt_f, phi_f, eta_f);
  return ((cGT.pt() == 0xFFFF) && (cGT.phi() == 0x1FFF) && (cGT.eta() == 0xF1D) && cGT.passNullBitsCheck());
}

/*
 * Check signed integer representation
 */
int passSignedIntegerRepr(void) {
  bool isValid = true;
  float pt_f = 2048;  // max value is 2047.97
  float phi_f = -1 * M_PI;
  float eta_f = -1.4841;
  float precision = 0.01;

  l1tp2::DigitizedClusterGT cGT = l1tp2::DigitizedClusterGT(isValid, pt_f, phi_f, eta_f);
  return (round_to(cGT.realEta(), precision) == round_to(eta_f, precision)) &&
         (round_to(cGT.realPhi(), precision) == round_to(phi_f, precision));
}

/*
 * Check signed integer representation
 */
int passSignedIntegerReprPositive(void) {
  bool isValid = true;
  float pt_f = 2048;  // max value is 2047.97
  float phi_f = -0.87266;
  float eta_f = 0.83;
  float precision = 0.01;

  l1tp2::DigitizedClusterGT cGT = l1tp2::DigitizedClusterGT(isValid, pt_f, phi_f, eta_f);
  return ((round_to(cGT.realEta(), 2) == round_to(eta_f, 2)) && (round_to(cGT.realPhi(), 2) == round_to(phi_f, 2)));
}

int main() {
  if (!passItemizedDigitizedConstructorTest()) {
    std::cout << "Fails itemized digitized constructor test!\n" << std::endl;
    return 1;
  }

  if (!passPositiveConstructorTest()) {
    std::cout << "Fails positive constructor test!\n" << std::endl;
    return 1;
  }

  if (!passNegativeConstructorTest()) {
    std::cout << "Fails negative constructor test!\n" << std::endl;
    return 1;
  }

  if (!passSignedIntegerRepr()) {
    std::cout << "Fails signed integer representation test!\n" << std::endl;
    return 1;
  }

  if (!passSignedIntegerReprPositive()) {
    std::cout << "Fails signed integer positive representation test!\n" << std::endl;
    return 1;
  }

  return 0;
}