% clear all
% clc
% %get the imaging strain
% fileNamecsv = strcat('strain');
% fileFormat  = '.csv';
% spotData = csvread( strcat(fileNamecsv,fileFormat) , 1 , 1 );
% Sdata    = spotData(:,1:3);
% %%%%%%%%%%%%%%%%%%%%%%%
% % --- define each spot
% BL = find(Sdata(:,1) == 0);
% BR = find(Sdata(:,1) == 1);
% TR = find(Sdata(:,1) == 2);
% TL = find(Sdata(:,1) == 3);
% 
% x(:,1) = Sdata(BL,2);
% x(:,2) = Sdata(BR,2);
% x(:,3) = Sdata(TR,2);
% x(:,4) = Sdata(TL,2);
% 
% y(:,1) = Sdata(BL,3);
% y(:,2) = Sdata(BR,3);
% y(:,3) = Sdata(TR,3);
% y(:,4) = Sdata(TL,3);
% 
% % --- x & y reference points
% xrefer = x(1,1:4);
% yrefer = y(1,1:4);
% 
% % --- Displacement
% u = bsxfun(@minus,x,xrefer); 
% v = bsxfun(@minus,y,yrefer); 
% 
% % --- Initialise
% Exx = zeros(length(x),1);
% Eyy = zeros(length(x),1);
% Exy = zeros(length(x),1);
% 
% % --- Shape function
% for i = 1:length(x)
%     xs = 1/4*(-x(i,1)+x(i,2)+x(i,3)-x(i,4));
%     ys = 1/4*(-y(i,1)+y(i,2)+y(i,3)-y(i,4));
%     xt = 1/4*(-x(i,1)-x(i,2)+x(i,3)+x(i,4));
%     yt = 1/4*(-y(i,1)-y(i,2)+y(i,3)+y(i,4));
%     
%     us = 1/4*(-u(i,1)+u(i,2)+u(i,3)-u(i,4));
%     ut = 1/4*(-u(i,1)-u(i,2)+u(i,3)+u(i,4));
%     vs = 1/4*(-v(i,1)+v(i,2)+v(i,3)-v(i,4));
%     vt = 1/4*(-v(i,1)-v(i,2)+v(i,3)+v(i,4));
%     
%     JJ = xs*yt-xt*ys;
%     
%     % --- Set up matrix
%     MA(1,1) = yt;
%     MA(1,2) = -ys;
%     MA(2,1) = -xt;
%     MA(2,2) = xs;
%     
%     MB1(1,1) = us;
%     MB1(2,1) = ut;
%     MB2(1,1) = vs;
%     MB2(2,1) = vt;
%     
%     MC1 = 1/JJ*MA*MB1;
%     MC2 = 1/JJ*MA*MB2;
%     
%     ux = MC1(1,1);
%     uy = MC1(2,1);
%     vx = MC2(1,1);
%     vy = MC2(2,1);
%     
%     %Calculate strain tensor
%     Exx(i) = abs(ux+1/2*(ux^2+vx^2));
%     Eyy(i) = abs(vy+1/2*(vy^2+uy^2));
%     Exx(i) = (Exx(i)+Eyy(i))/2-0.075-0.138-0.02;
%     Eyy(i) = Exx(i);
%     Exy(i) = abs(1/2*(uy+vx+ux*uy+vx*vy));
% end
%%%%%%%%%%%%%%%%%%%%
side = 5/1000;
thick = 0.8/1000;
Exx = 0:0.001:1;
Eyy = 0:0.001:1;
Tx = f7(Exx);
Ty = f8(Eyy);
Fx=Tx*side*thick;
Fy=Ty*side*thick;
% for i = 3
%     Fx = smooth(Fx);
%     Fy = smooth(Fy);
%     Exx = smooth(Exx);
%     
% end
%%%%%%%%%%%%%%%%%%%%
% FxName = strcat('Mainx');
% FyName = strcat('Mainy');
% biaxialFormat = '.txt';
% Fx = textread(strcat(FxName,biaxialFormat));
% Fx = -Fx(:,2);
% Fx = Fx - min(Fx);
% Fy = textread(strcat(FyName,biaxialFormat));
% Fy = -Fy(:,2);
% Fy = Fy - min(Fy);
%%%%%%%%%%%%%%%%%%%%%
FxNameMPM = strcat('DSMx');
FyNameMPM = strcat('DSMy');
biaxialFormat = '.txt';
FxMPM = textread(strcat(FxNameMPM,biaxialFormat));
FxMPM = -FxMPM(:,2);
FxMPM = FxMPM- min(FxMPM);
FyMPM = textread(strcat(FyNameMPM,biaxialFormat));
FyMPM = -FyMPM(:,2);
FyMPM = FyMPM- min(FyMPM);
figure; 
plot(FyMPM)
   [x_del y_del] = ginput();
    x_del         = round(x_del);
    
for i = 1:length(y_del)
    for j = 1:length(Fy)
        tmp(j) = abs(Fy(j)-y_del(i));   
    end
     [idy idx] = min(tmp);
        straincaly(i) = Eyy(idx)
end

figure(2); 
plot(FxMPM)
   [x_del y_del] = ginput();
    x_del         = round(x_del);
    
for i = 1:length(y_del)
    for j = 1:length(Fx)
        tmp(j) = abs(Fx(j)-y_del(i));   
    end
     [idy idx] = min(tmp);
        straincalx(i) = Exx(idx);
end

straincal = (straincalx+straincaly)/2
%save MPMlumen straincal    
    
    
    
    
    
    
    
    
    
    