% -------------------------------------------------------------------------
%                    Manual Strain + biaxial processing
% -------------------------------------------------------------------------
% Loads .csv file from spot tracker using 4 dots (counter clockwise from
% bottom left). Calculates biaxial strain and plots againts stress from
% load cell.

% Choose sample number below. Arrays give measured dimensions of each
% samples with the offset for each curve to match the stress/strain peaks.
% comment out the data type that is not being used, i.e. FY, FA, RB
clear all, close all

% --- Uncomment for plotting at the end
% export_fig biaxial_FA04.jpg -transparent

% --- Turn biaxial data plotting on or off if data exists
manual  = 1;
biaxial = 1;
sample  = 1;

% samples 4,5,7 bad strain measurements (double peaks)

% --- Specify file names and format
fileNamecsv = strcat('strain');
fileFormat  = '.csv';

% --- Sample Dimensions and peak matches RB
%         RB01 RB02 RB03 RB04 RB05 RB06 RB07 RB08 RB09
% x_d    = [7.0, 6.0, 6.0, 6.0, 4.0, 5.0, 5.0, 5.0, 5.0]/1000; 
% y_d    = [7.0, 6.0, 5.5, 6.0, 4.0, 5.0, 5.0, 5.0, 5.0]/1000;
% z_d    = [1.2, 1.1, 1.2, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0]/1000;
% E_peak = [92,  159, 222, 70,  86,  83,  39,  148, 153];
% S_peak = [126, 168, 224, 122, 58,  56,  67,  150, 178];

% --- Sample Dimensions and peak matches FY
%      FY01 FY02
% x_d = [5.5, 5.5]/1000; 
% y_d = [5.5, 5.5]/1000;
% z_d = [1.2, 1.2]/1000;
% 
% E_peak = [0, 0];
% S_peak = [0, 0];

% --- Sample Dimensions and peak matches FA
%      FA01 FA02 FA03 FA04
x_d = [5]/1000; 
y_d = [5]/1000;
z_d = [0.8]/1000;

E_peak = [62];
S_peak = [60];

% --- define first peak coordinate manually to match curves
peak = abs(S_peak(sample)) - abs(E_peak(sample));

% --- Load tracking data
spotData = csvread( strcat(fileNamecsv,fileFormat) , 1 , 1 );
Sdata    = spotData(:,1:3);

% --- Load biaxial data if biaxial = 1
if biaxial == 1
    FxName = strcat('mainx');
    FyName = strcat('mainy');
    biaxialFormat = '.txt';
    Fx = textread(strcat(FxName,biaxialFormat));
    Tx = - Fx(1:1:length(Fx),2) / (x_d(sample) * z_d(sample));
    
    Fy = textread(strcat(FyName,biaxialFormat));
    Ty = - Fy(1:1:length(Fy),2) / (y_d(sample) * z_d(sample));
    
    if manual == 1
        % --- Match peaks
        Tx(1:peak) = [];
        Ty(1:peak) = [];
    end

end

% --- define each spot
BL = find(Sdata(:,1) == 0);
BR = find(Sdata(:,1) == 1);
TR = find(Sdata(:,1) == 2);
TL = find(Sdata(:,1) == 3);

x(:,1) = Sdata(BL,2);
x(:,2) = Sdata(BR,2);
x(:,3) = Sdata(TR,2);
x(:,4) = Sdata(TL,2);

y(:,1) = Sdata(BL,3);
y(:,2) = Sdata(BR,3);
y(:,3) = Sdata(TR,3);
y(:,4) = Sdata(TL,3);

% --- x & y reference points
xrefer = x(1,1:4);
yrefer = y(1,1:4);

% --- Displacement
u = bsxfun(@minus,x,xrefer); 
v = bsxfun(@minus,y,yrefer); 

% --- Initialise
Exx = zeros(length(x),1);
Eyy = zeros(length(x),1);
Exy = zeros(length(x),1);

% --- Shape function
for i = 1:length(x)
    xs = 1/4*(-x(i,1)+x(i,2)+x(i,3)-x(i,4));
    ys = 1/4*(-y(i,1)+y(i,2)+y(i,3)-y(i,4));
    xt = 1/4*(-x(i,1)-x(i,2)+x(i,3)+x(i,4));
    yt = 1/4*(-y(i,1)-y(i,2)+y(i,3)+y(i,4));
    
    us = 1/4*(-u(i,1)+u(i,2)+u(i,3)-u(i,4));
    ut = 1/4*(-u(i,1)-u(i,2)+u(i,3)+u(i,4));
    vs = 1/4*(-v(i,1)+v(i,2)+v(i,3)-v(i,4));
    vt = 1/4*(-v(i,1)-v(i,2)+v(i,3)+v(i,4));
    
    JJ = xs*yt-xt*ys;
    
    % --- Set up matrix
    MA(1,1) = yt;
    MA(1,2) = -ys;
    MA(2,1) = -xt;
    MA(2,2) = xs;
    
    MB1(1,1) = us;
    MB1(2,1) = ut;
    MB2(1,1) = vs;
    MB2(2,1) = vt;
    
    MC1 = 1/JJ*MA*MB1;
    MC2 = 1/JJ*MA*MB2;
    
    ux = MC1(1,1);
    uy = MC1(2,1);
    vx = MC2(1,1);
    vy = MC2(2,1);
    
    %Calculate strain tensor
    Exx(i) = abs(ux+1/2*(ux^2+vx^2));
    Eyy(i) = abs(vy+1/2*(vy^2+uy^2));
    Exx(i) = (Exx(i)+Eyy(i))/2;
    Eyy(i) = Exx(i);
    Exy(i) = abs(1/2*(uy+vx+ux*uy+vx*vy));
end

if manual == 1
    % --- Manually pick regions to delete (noise, etc) (if manual = 1)
    % --- Select first and last points of each test in sequence
    
    figure; plot(Exx)
    [x_del y_del] = ginput();
    x_del         = round(x_del);
    
    figure;
    for i = 1:length(x_del)/2
        
        % --- Define each loading cycle individually
        Sx{i} = Tx(x_del((i*2)-1):x_del(i*2))';
        Sx{i} = Sx{i} - Tx(x_del((i*2)-1));
        

        Sy{i} = Ty(x_del((i*2)-1):x_del(i*2))';
        Sy{i} = Sy{i} - Ty(x_del((i*2)-1));

        Ex{i} = Exx(x_del((i*2)-1):x_del(i*2))';
        Ey{i} = Eyy(x_del((i*2)-1):x_del(i*2))';
        %zero the first data
        Ex{i} = Ex{i} - Ex{i}(1);
        Ey{i} = Ey{i} - Ey{i}(1);

        plot(Ex{i},Sx{i},'r*','lineWidth',1.5), hold on
        plot(Ey{i},Sy{i},'b*','lineWidth',1.5)
        xlabel('Strain','fontSize',14), ylabel('Stress (Pa)','fontSize',16)
        legend('Circumferential','Longitudinal','Location','NW')
        
         end
else    
    % --- Plot stress vs. strain if biaxial = 1
    if biaxial == 1
        Tx = Tx - min(Tx);
        Ty = Ty - min(Ty);
        
        figure;
        plot(Exx,Tx,'r','lineWidth',1.5), hold on
        plot(Eyy,Ty,'b','lineWidth',1.5)
        xlabel('Strain','fontSize',14), ylabel('Stress (Pa)','fontSize',16)
        legend('Circumferential','Longitudinal','Location','NW')
    else
        % --- Strain figures
        figure;
        plot(Exx,'r.','lineWidth',1.5), hold on
        plot(Eyy,'b.','lineWidth',1.5)
        xlabel('frames','fontSize',14), ylabel('\epsilon','fontSize',16)
        legend('Ex','Ey','Location','NW')
    end
end

% --- Plots for checking peaks for matching
figure;
subplot(1,2,1)
plot(Ty)
title('Stress')
subplot(1,2,2)
plot(Eyy)
title('Strain')
% --- Plots for checking peaks for matching
figure;
subplot(1,2,1)
plot(Tx)
title('Stress')
subplot(1,2,2)
plot(Exx)
title('Strain')
