clc
clear
load('INFO_FE.mat');
load('dataP4.mat');

neig = 21;
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
e = 2.71;

displacement = zeros(size(MODES, 1), timeStep); %displacement en cada node, time quan passa 
t = 0;
cont = 1;
Amplitud = zeros(neig, size(MODES, 1));
for j = 1:1:timeStep
    dTime = 0;
    for iMode = 1:1:neig
        wi = 2*pi*FREQ(iMode)*sqrt(1-dampingRatio^2);
        qo = MODES(:,iMode)'*Mll*dll;
        newDTime = (MODES(iMode)*(e^(-dampingRatio*2*pi*FREQ(iMode)*t)*(qo*cos(wi*t)+sin(wi*t)*(dampingRatio*qo)/sqrt(1-dampingRatio^2))))';
        

        if j == 1
            Amplitud(iMode, cont) = MODES(:,iMode)'*Mll*newDTime; % 
            cont = cont +1;
        end
        

        dTime = dTime + newDTime;
    end
    displacement(:, j) = dTime; 
    t = t + time/timeStep;
end

plot(FREQ, Amplitud(: , 2))


