%plots for the paper
figure
subplot(1,2,1)
t= 0:.001:2;
a = 1;
f = 7;
p = pi/2;
y = a*sin(2*pi*f*t + p);
plot(t,y);
xlabel('Time (s)'), %ylabel('Amplitude')
title('fq =7')
axis([0 2 -a*2 a*2])

subplot(1,2,2)
t= 0:.001:2;
a = 1;
f = 2;
p = 0;
y = a*sin(2*pi*f*t + p);
plot(t,y);
xlabel('Time (s)'), %ylabel('Amplitude')
title('fq =2')
axis([0 2 -a*2 a*2])