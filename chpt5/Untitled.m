itload('fft_file.it')  % load the exported it++ file
figure;
x=1:512;
y=1:601;
plot(x,freq_signal,'b','LineWidth',2);
grid;