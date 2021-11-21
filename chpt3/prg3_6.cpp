#include <itpp/itcomm.h>
  using namespace itpp;

  int main() {
// Initiate the BSC with cross-over probability 0.1
    BSC bsc(0.1);

    bvec transmitted_bits = randb(100);
    bvec received_bits = bsc(transmitted_bits);
  }
