
global G_const G_limits betatarg RBtype linStressComputed linear_prob
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl DinamcVar Method
global serieLimitState Stochast_moment GradOutput_on
linStressComputed=0;
GradOutput_on=1;

SpaceTruss_25
 
arq.x0=[1 1 1 1 1 1]*2.3;

% %Article Min mean(c)
%  arq.x0=[0.0500	0.0500	5.74	1.718	1.054	5.574];
% %Article Min std(c)
 arq.x0=[0.147	0.672	3.465	0.566	0.822	8.048];


% 
% arq.x0=[.1 .1 1 1 2 6];
% 
% % 
  
  X_ops=[0.2420 0.3087  3.5243  0.5693  1.0858  7.6058
        0.0718	0.1360	4.0281	0.5548	1.1669	7.2030
        0.0500	0.0500	4.3690	0.5406	1.2682	6.7505
        0.0500	0.0500	4.7735	0.5271	1.2987	6.3515
        0.0500	0.0500	4.8287	0.6522	1.2652	6.3279
        0.0500	0.0500	4.9631	0.8469	1.1292	6.4801
        0.0500	0.0500	5.0234	1.0567	1.0938	6.4101
        0.0500	0.0500	5.1304	1.3117	1.0563	6.2830
        0.0500	0.0500	5.4264	1.4534	1.0471	5.9915
        0.0500	0.0500	5.6769	1.6588	1.0457	5.6795];
%       0.0500  0.0500  4.4830  1.6792  1.4135  5.6249

% x =
% 
%     0.0500    0.0500    4.7882    1.8365    1.1051    6.1295
%     0.0500    0.0500    4.4831    1.6790    1.4135    5.6249
%     0.0500    0.0500    4.4826    1.8637    1.0658    6.4653
%     0.0500    0.0500    4.5092    1.8288    1.1338    6.2792
%     0.0500    0.0500    4.5474    1.8299    1.1297    6.2591
%     0.0500    0.0500    4.5862    1.8311    1.1256    6.2386
%     0.0500    0.0500    4.6255    1.8322    1.1215    6.2177
%     0.0500    0.0500    4.6654    1.8333    1.1174    6.1963
%     0.0500    0.0500    4.7058    1.8344    1.1133    6.1744
%     0.0500    0.0500    4.7467    1.8354    1.1092    6.1522
% 
% 
%   X_ops=[0.0500    0.0500    5.6752    1.6600    1.0455    5.6806
%     0.0500    0.0500    5.2865    1.3657    1.0568    6.1263
%     0.0500    0.0500    4.9685    1.0556    1.0869    6.4733
%     0.0500    0.0500    4.7074    0.6777    1.1902    6.6144
%     0.0500    0.0729    4.4427    0.5442    1.2512    6.7132
%     0.0500    0.1566    4.1089    0.5399    1.2178    6.9915
%     0.4351    1.1435    2.7577    0.3445    1.1073    7.4633
%     1.0107    2.0972    1.2531    0.0795    0.9196    8.3123
%     1.0176    2.0439    0.8496    0.0500    0.7589    9.1393
%     0.8945    1.8315    0.5940    0.0500    0.6057    9.9847];

  %arq.x0=X_ops;
%   arq.x0=X_ops(1,:);

X_ops = arq.x0;


%RV_consider = [6,11,13,14];
%DinamcVar = 3;Method='PC';

% 
% [x,f/1000]
% 
% ans =
% 
%   0.0500  0.0500  5.6769  1.6588  1.0457  5.6795  5.3257  0.3752
%   0.0635  0.8240  2.7811  0.0515  1.2121  7.6868  7.0011  0.2527(t25 s.27)

%new0.0728  0.3072  3.5624  0.1907  1.2889  7.2319  6.5440  0.2642

%   0.3544  1.0157  2.9547  0.0500  1.1908  7.3809  6.6462  0.2540(t23 s.27)
%   0.3068  0.9230  3.3664  0.0500  1.3464  6.7264  6.3631  0.2605
%   0.3048  0.9473  3.6826  0.0500  1.4616  6.1372  6.1312  0.2708(t20 s.2830)
%   0.0500  0.1844  4.4947  0.2951  1.4243  6.2326  5.9377  0.2839
%   0.0500  0.0500  4.6861  0.5221  1.2814  6.4711  5.7510  0.2974
%   0.0500  0.0500  4.7957  0.7857  1.1547  6.5794  5.5835  0.3124
%   0.0500  0.0500  4.9996  1.0804  1.0837  6.4432  5.4472  0.3296
%   0.0500  0.0500  5.3009  1.3761  1.0557  6.1119  5.3574  0.3503
%arq.x0=[0.0500  0.0500  5.6769  1.6588  1.0457  5.6795];

%arq.x0=[.6132    1.3765    2.0653    0.1500    1.0461    8.0964];

%  arq.x0=[1 1 2 1 1 3];
  %arq.x0=[1 1 5 1 1 5];
 %arq.x0=[.05 .05 5.74 1.72 1.05 5.57];
 %SD opt [6,9:14][ 0.3788    0.1717    3.8639    0.0526    1.7117    6.0564]
 %[0.0635    0.8240    2.7811    0.0515    1.2121    7.6868]
 
%arq.x_paper= [.147, .672, 3.465, 0.566, 0.822, 8.048];
%arq.x0=[.147, .672, 3.465, 0.566, 0.822, 8.048];

 
%N/mm^2 = MPa = 1MN/(1e2cm)^2 = 1MN/1e4cm^2 = 1e2N/cm^2
 %2e5MPa = 2e7 N/cm^2

Stochast_moment = 1;
% DinamcVar = 5000;Method='MC';
% DinamcVar = 2;Method='MC';
DinamcVar = 2;Method='PC';

%DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';


%k_p = 2000;%500
k_p = 1;
RBtype=1;%POD
%RBtype=2;%RBM
 RBtype=0;%FEM
linear_prob=1;
%RBtype=0;
 
 %% Constraint
global ResType n_sigma sign_const

%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy

arq.tpres=[0 0 0 0 0 0];% Reliability
 
 iconst = 3;%stress
 arq.tpres(iconst)=1;
 ResType{iconst}='probability';
 n_sigma(iconst)=3;
 sign_const(iconst)=1;
 arq.cm = -15000;
 arq.em =  5000;
 

iconst = 1;%vol
arq.tpres(iconst)=1;
ResType{iconst}='deterministic';
n_sigma(iconst)=0;
sign_const(iconst)=1;
 arq.vol0 = 750;

%% Objective
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd

n_tot_obj = 5;
global mo  obj_type i_obj
mo = zeros(1,n_tot_obj*2);
obj_type=cell(1,n_tot_obj*2);

main_object=4;% energy
arq.tobj=0;%Multi_Obj
mo([main_object,n_tot_obj+main_object])=1;
 obj_type{main_object}='single';
 obj_type{n_tot_obj+main_object}='single';
 %i_obj = 661; %specific
 arq.tobj=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarDist (6,nrv) = [ distribuction; mean; std; par3; par4;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

Stochast_on=1;
nRandVar = 14;
var_tpe(1:6) = 2; % Sensibility variable E
var_tpe(7:8) = 3; % Sensibility variable Px
var_tpe(9:14) = 1; % Sensibility variable x

E=arq.E;dE=E/5;% 5,0%


%Random Variables
% 1 - E all design regions (3)
% 2 - Load (1)
% RandVarProp(1,:)=[2 2 2 3]; % type var
% RandVarProp(2,:)=[1 1 1 5]; % which type data
% rv_data_comp = {1,2,3,0};
% 
% RandVarDist(1,:)=[2 2 2 1]; % distribuction
% RandVarDist(2,:)=[E E E 1]; % Mu par1
% RandVarDist(3,:)=[dE dE dE .25]; % std par2


RandVarProp(1,:)=var_tpe; % type var
RandVarProp(2,1: 6)=1; % which type data
RandVarProp(2,7: 8)=4; % which type data 
RandVarProp(2,9:14)=5; % which type data

%(DoF node i direction x = i*3-2)
rv_data_comp = {1,2,3,4,5,6,3*3-2,6*3-2,1,2,3,4,5,6};

RandVarDist(1,:)=1*ones(1,nRandVar); % distribuction
RandVarDist(2,1: 6)=1e7; % Mu par1
RandVarDist(2,7: 8)=500; % Mu par1
RandVarDist(2,9:14)=1; % Mu par1

RandVarDist(3,1: 5)=2e5; % std par2
RandVarDist(3,   6)=1.5e6; % std par2
RandVarDist(3,7: 8)=50; % std par2
RandVarDist(3,9:14)=0.05; % std par2

%%%%%%%%%%%%%%%%%%%%%%%%
% R. V. Elimination!! 
% R. V. Considered:
RV_consider = [5,6,11,13,14];
% RV_consider = [1:14];
nRandVar = length(RV_consider);
rv_data_comp = {rv_data_comp{RV_consider}};
RandVarDist=RandVarDist(:,RV_consider);
RandVarProp=RandVarProp(:,RV_consider);
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit State (G=0)
G_fun = 'StructLimitState';
S_fun = 'Stoch_Analysis';
%G_const = 1 VOLUME, 2 FLAMBAGEM GLOBAL, 3 FLAMBAGEM LOCAL, 4 TENSAO, 
% 5 DESLOCAMENTO, 7 DESLOCAMENTO ESPECIFICO, 8 TENSAO ESPECIFICA
%G_const = [3,4];
%G_limits = [1,10e3];
% G_const = 4;
% G_limits = 65000;

 G_const = 5;
 G_limits = 1.2;
 
serieLimitState = [];

betatarg = 3.;
%Mu = [1,1,1]*.5;
%Cu = diag([1 2 3]/50);
%ft=[2,2,2];


nrv=size(RandVarProp,2);
Correl=eye(nrv);


arq=MatProp(arq);
%SIGu=arq.MatCurv(end)*.5;