function [stimuli] = CMR(CF,FM,FBs)

T = 0.5;
T2 = 0.05;
sr = 100e3;
dt = 1/sr;
t = 0:dt:T-dt;
t2 = 0:dt:T2-dt;
m = 1;




OFM = sin(2*pi*CF.*t).*(1+m*(cos(2*pi*FM.*t+pi)));

S = sin(2*pi*CF/2.*t2);
Gate = sin(2*pi*10.*t2).^2;
S = 2*(S.*Gate);
Spad = [zeros(1,3.5*length(S)) S zeros(1,length(S)) S zeros(1,length(S)) S zeros(1,1.5*length(S))];

for i = 1:length(FBs)
    CMFB(i,:)  = sin(2*pi*FBs(i).*t).*(1+m*(cos(2*pi*FM.*t+pi)));
    CDFB(i,:)  = sin(2*pi*FBs(i).*t).*(1+m*(cos(2*pi*FM.*t)));
end


foo;

