clc
clear
load('INFO_FE.mat');
load('dataP4.mat');

neig = 25;
time = 105.23; % 1/0.3801*40
timeStep = 500;
dampingRatio = 0.01;

nnode = size(COOR,1);
ndim = size(COOR,2);
DOFl = 1:nnode*ndim;
DOFl(DOFr) = [];

Mll = M(DOFl,DOFl);
Kll = K(DOFl,DOFl);
dll = d(DOFl);
[MODES, FREQ] = UndampedFREQ(Mll,Kll,neig);
ee = 2.71;

displacement = zeros(size(MODES, 1), timeStep); %displacement en cada node, time quan passa
t = 0;
cont = 1;
Amplitud = zeros(neig, 1);
for j = 1:1:timeStep
    dTime = 0;
        for iMode = 1:1:neig
            wi = FREQ(iMode)*sqrt(1-dampingRatio^2);
            qo = MODES(:,iMode)'*Mll*dll;
            q = ee^(-dampingRatio*FREQ(iMode)*t)*(qo*cos(wi*t)+sin(wi*t)*(dampingRatio*qo)/sqrt(1-dampingRatio^2));
            phy = MODES(:, iMode);
            dis = phy*q;
            dTime = dTime + dis;
            if j == 1
                Amplitud(iMode) = abs(qo); %calcula l'amplitud
            end
        end
    displacement(:, j) = dTime; % obtain eq 95
    t = t + time/timeStep;
end

bar(Amplitud(:))
xlabel('Modes')
ylabel('Amplitude')
