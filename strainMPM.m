
Fxx = Fx(:,2);
Fyy = Fy(:,2);
Fxx = Fxx(1:2:length(Fxx)-1)
Fyy = Fyy(1:2:length(Fyy)-1)
Exx = Exx(1:length(Exx)-1)
Eyy = Eyy(1:length(Eyy)-1)
for i = 1: 2
Fxx = smooth(Fxx);
Fyy = smooth(Fyy);
end
figure
plot(Exx,Fxx)
figure(2)
plot(Eyy,Fyy)
