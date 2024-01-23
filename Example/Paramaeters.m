%% Configure Simulation
%Define the nominal plant model.

C = 0.5;
L = 0.5;
T = 0.5;
G = tf(C,[T 1],'Outputdelay',L);

%Generate the conventional PID controller parameters using pidtune.

pidController = pidtune(G,'pidf'); 

%In this example, the reference (r) is a step signal and t_r = 0, 
% which results in Ce=1 as follows: 
Ce = 1;
%To configure the simulation, use the following nominal controller parameters.

tauC = 0.2;

Cd = min(T,L/2)*Ce;
C0 = 1/(C*Ce*(tauC+L/2));
C1 = max(T,L/2)*C0; 


%% Coonctruct Fuzzy inference systems 

%Create a type-1 FIS using sugfis.
fis1 = sugfis;
%Add input variables to the FIS.
fis1 = addInput(fis1,[-1 1],'Name','E');
fis1 = addInput(fis1,[-1 1],'Name','delE');

%Add three uniformly distributed overlapping triangular membership functions, 
% (MFs) to each input. The MF names stand for negative (N), zero (Z), 
% and positive (P).
fis1 = addMF(fis1,'E','trimf',[-2 -1 0],'Name','N');
fis1 = addMF(fis1,'E','trimf',[-1 0 1],'Name','Z');
fis1 = addMF(fis1,'E','trimf',[0 1 2],'Name','P');
fis1 = addMF(fis1,'delE','trimf',[-2 -1 0],'Name','N');
fis1 = addMF(fis1,'delE','trimf',[-1 0 1],'Name','Z');
fis1 = addMF(fis1,'delE','trimf',[0 1 2],'Name','P');

%Plot the input membership functions.
% figure
% subplot(1,2,1)
% plotmf(fis1,'input',1)
% title('Input 1')
% subplot(1,2,2)
% plotmf(fis1,'input',2)
% title('Input 2')

%Add the output variable to the FIS.
fis1 = addOutput(fis1,[-1 1],'Name','U');
%Add uniformly distributed constant functions to the output, 
% The MF names stand for negative big (NB), negative medium (NM), zero (Z),
% positive medium (PM), and positive big (PB).

fis1 = addMF(fis1,'U','constant',-1,'Name','NB');
fis1 = addMF(fis1,'U','constant',-0.5,'Name','NM');
fis1 = addMF(fis1,'U','constant',0,'Name','Z');
fis1 = addMF(fis1,'U','constant',0.5,'Name','PM');
fis1 = addMF(fis1,'U','constant',1,'Name','PB');

%Add rules to the FIS. These rules create a proportional control surface.
rules = [...
    "E==N & delE==N => U=NB"; ...
    "E==Z & delE==N => U=NM"; ...
    "E==P & delE==N => U=Z"; ...
    "E==N & delE==Z => U=NM"; ...
    "E==Z & delE==Z => U=Z"; ...
    "E==P & delE==Z => U=PM"; ...
    "E==N & delE==P => U=Z"; ...
    "E==Z & delE==P => U=PM"; ...
    "E==P & delE==P => U=PB" ...
    ];
fis1 = addRule(fis1,rules);

%Convert the type-1 FIS, fis1, to a type-2 FIS.
fis2 = convertToType2(fis1);
scale = [0.2 0.9 0.2;0.3 0.9 0.3];
for i = 1:length(fis2.Inputs)
    for j = 1:length(fis2.Inputs(i).MembershipFunctions)
        fis2.Inputs(i).MembershipFunctions(j).LowerLag = 0;
        fis2.Inputs(i).MembershipFunctions(j).LowerScale = scale(i,j);
    end
end



%Plot the type-2 input membership functions.

% figure
% subplot(1,2,1)
% plotmf(fis2,'input',1)
% title('Input 1')
% subplot(1,2,2)
% plotmf(fis2,'input',2)
% title('Input 2')
