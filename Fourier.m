
rij=load('rijThetaFFT_tersoff-iter100050')
r=rij.rijFFT;
size(r)

f=1/1.018e-14*(0:256)/512;

y12=fft(r(:,1),512);
pyy12=y12.*conj(y12)/512;
subplot(3,3,1);
plot(f(2:257),pyy12(2:257));

y13=fft(r(:,2),512);
pyy13=y13.*conj(y13)/512;
subplot(3,3,2);
plot(f(2:257),pyy13(2:257));

y14=fft(r(:,3),512);
pyy14=y14.*conj(y14)/512;
subplot(3,3,3);
plot(f(2:257),pyy14(2:257));

y15=fft(r(:,4),512);
pyy15=y15.*conj(y15)/512;
subplot(3,3,4);
plot(f(2:257),pyy15(2:257));


y312=fft(r(:,5),512);
pyy312=y312.*conj(y312)/512;
subplot(3,3,5);
plot(f(2:257),pyy312(2:257));

y412=fft(r(:,6),512);
pyy412=y412.*conj(y412)/512;
subplot(3,3,6);
plot(f(2:257),pyy412(2:257));

y512=fft(r(:,7),512);
pyy512=y512.*conj(y512)/512;
subplot(3,3,7);
plot(f(2:257),pyy512(2:257));



z12= trapz(f(2:257),pyy12(2:257))
z13= trapz(f(2:257),pyy13(2:257))
z14= trapz(f(2:257),pyy14(2:257))
z15= trapz(f(2:257),pyy15(2:257))
z312=trapz(f(2:257),pyy312(2:257))
z412=trapz(f(2:257),pyy412(2:257))
z512=trapz(f(2:257),pyy512(2:257))


pyy12_NORMALIZED= pyy12./z12;
pyy13_NORMALIZED= pyy13./z13;
pyy14_NORMALIZED= pyy14./z14;
pyy15_NORMALIZED= pyy15./z15;
pyy312_NORMALIZED= pyy312./z312;
pyy412_NORMALIZED= pyy412./z412;
pyy512_NORMALIZED= pyy512./z512;


pyy_ADD = pyy12 + pyy13 + pyy14 + pyy15 + pyy312 + pyy412 + pyy512;
