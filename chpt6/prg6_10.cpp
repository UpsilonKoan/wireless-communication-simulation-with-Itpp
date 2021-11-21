#include <iomanip>
#include <itpp/itcomm.h>
using namespace std;
using namespace itpp;

int main()
{  
  ND_UPAM chan; // multidimensional channel with uniform PAM
  chan.set_M(3, 4); // 3-dimensional matrix channel, 4-PAM per dimension
  cout << chan << endl;
  bvec b = randb(3*2); // 3*2 bits in total
  vec x = chan.modulate_bits(b);
  mat H = randn(4,3); // 4 x 3 real matrix channel3个发射天线 4个接收天线
  double sigma2 = 0.01; // noise variance per real dimension
  vec y = H*x + sqrt(sigma2)*randn(4); // transmit vector x
  QLLRvec llr; // log-likelihood ratios
  QLLRvec llr_ap = zeros_i(3*2);  // a priori equiprobable bits
  chan.demodulate_soft_bits(y, H, sigma2, llr_ap, llr);
  cout << "True bits:" << b << endl;
  cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
}
