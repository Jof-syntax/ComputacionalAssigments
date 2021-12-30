clc
clear
load('INFO_FE.mat');        % these loads contain some useful information
load('dataP4.mat');         % such as the COOR and DORr arrays or
                            % the stiffness and mass matrices

neig = 25;                  % number of modes
timeStep = 500;             % number of time steps, imposed by section 8 of the assigment
dampingRatio = 0.01;        % dumping ratio (1%)

nnode = size(COOR,1);
ndim = size(COOR,2);
DOFl = 1:nnode*ndim;
DOFl(DOFr) = [];

Mll = M(DOFl,DOFl);
Kll = K(DOFl,DOFl);
dll = d(DOFl);
[MODES, FREQ] = UndampedFREQ(Mll,Kll,neig);         %using this function
                                                    %we obtain the natural
                                                    %frequencies and modes
time = 40*2*pi/FREQ(1); % Computes the maximum time as 40 time the period
ee = 2.71; % Mathematical constant

DISP = zeros(size(d,1), timeStep); % Generates the displacement matrix, where the rows will be the displacements and the columns each time.
t = zeros(timeStep, 1); % Initial time 
for j = 1:1:timeStep-1 % Obtains the 't' matrix
        t(j+1) = t(j) + time/timeStep;  
end
Amplitud = zeros(neig, 1); % Generates an amplitude matrix
for j = 1:1:timeStep % Obtains the 'DISP' matrix and compute the amplitude
    dTime = 0;
        for iMode = 1:1:neig
            wi = FREQ(iMode)*sqrt(1-dampingRatio^2); % Natural 
            qo = MODES(:,iMode)'*Mll*dll; % Initial amplitude
            q = ee^(-dampingRatio*FREQ(iMode)*t(j))*(qo*cos(wi*t(j))+sin(wi*t(j))*(dampingRatio*qo)/sqrt(1-dampingRatio^2)); %Modal coordinates
            phy = MODES(:, iMode); % Associated Mode
            dis = phy*q;
            dTime = dTime + dis; 
            if j == 1
                Amplitud(iMode) = abs(qo); % Computes the amplitude for t = 0.
            end
        end
    DISP([1:5760], j) = dTime; % Eq 95 of the theory pdf
end


 bar(Amplitud(:)) % Plot of the amplitud for each mode
 xlabel('Modes')
 ylabel('Amplitude')

GidPostProcessDynamic( COOR, CN, TypeElement, DISP , 'dataP4.mat', posgp' , 'MyFirstMesh3D.msh', t);
