%Get the MPM force data
FxName = strcat('MPM lumen x');
FyName = strcat('MPM lumen y');
biaxialFormat = '.txt';
MPMfx = textread(strcat(FxName,biaxialFormat));
MPMfy = textread(strcat(FyName,biaxialFormat));
figure2
plot(MPMfx)
figure3
plot(MPMfy)