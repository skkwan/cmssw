#include "DataFormats/L1TCorrelator/interface/TkHT.h"

using namespace l1t;

TkHT::TkHT() {}

TkHT::TkHT(double etTotal,
                   const edm::RefProd<TkJetCollection>& jetCollRef,
                   const edm::Ref<VertexCollection>& avtxRef,
                   int bx)
    : etTot_(etTotal), jetCollectionRef_(jetCollRef), vtxRef_(avtxRef), bx_(bx) {
  if (vtxRef_.isNonnull()) {
    setVtx(vtxRef()->z0());
  }
}
