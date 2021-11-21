#include <itpp/itcomm.h>
#include <iomanip>
using namespace std;
using namespace itpp;

int main()
{  
  ND_UQAM chan; // multidimensional channel with uniform QAM
  chan.set_M(3, 4); // 3-dimensional matrix channel, QAM per dimension
  cout << chan << endl;
  bvec b = randb(3*2); // 3*2 bits in total
  cvec x = chan.modulate_bits(b);
  cmat H = randn_c(4,3); // 4 x 3 real matrix channel
  double sigma2 = 0.01; // noise variance per real dimension
  cvec y = H*x + sqrt(sigma2)*randn_c(4); // transmit vector x
  QLLRvec llr; // log-likelihood ratios
  QLLRvec llr_ap = zeros_i(3*2);  // a priori equiprobable bits
  chan.demodulate_soft_bits(y, H, sigma2, llr_ap, llr);
  cout << "True bits:" << b << endl;
  cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
 }
