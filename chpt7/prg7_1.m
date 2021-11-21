figure(1); clf;
itload('ofdm_result_file.it');
h1 = semilogy(EbN0dB,ber,'*-r'); hold on
ebn0db = 0:.1:10;
ebn0 = 10.^(ebn0db/10);
Pb = 0.5 * erfc(sqrt(ebn0));
h2 = semilogy(ebn0db,Pb);
set(gca,'fontname','times','fontsize',16);
xlabel('{\it E_b} / {\it N}_0 [dB]');
ylabel('BER')
title('OFDM on an AWGN Channel');
legend([h1 h2],'Simulation','Theory');
grid on;
