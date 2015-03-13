
global G_const G_limits betatarg RBtype linStressComputed linear_prob
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl DinamcVar Method
global serieLimitState Stochast_moment var_type GradOutput_on
linStressComputed=0;
GradOutput_on=0;
Stochast_on=0;

% Created by Renato S. Motta (03/2014)
global bindir HYPLASexe
global NO_IG HYPLAS_flag hyplasprojname
NO_IG=1;
HYPLAS_flag=1;
bindir = 'bin';
HYPLASexe = 'hyplas90a.exe';

hyplasprojname='13_6_2';
 
arq.x0=[4,4];

arq.xmin=1;
arq.xmax=6;

nvar = length(arq.x0);
arq.arex = 1:nvar;
%N/mm^2 = MPa = 1MN/(1e2cm)^2 = 1MN/1e4cm^2 = 1e2N/cm^2
 %2e5MPa = 2e7 N/cm^2

Stochast_moment = 0;
DinamcVar = 5000;Method='MC';
%DinamcVar = 3;Method='PC';

%DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';


%k_p = 2000;%500
k_p = 1;
RBtype=1;%POD1
%RBtype=2;%RBM2
%RBtype=0;%FEM0
linear_prob=0;
%RBtype=0;

%% Objective
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};

        %arq.tobj=2;%    1---VOLUME  2----ENERGIA   3--Desl.Tot
        %4--Desl.Maxim   5--Desl.Esp   6--Tes.Max
        %arq.tpres=1;%    1---PESO    2----TENSAO     3---Tens e Desloc
        %4---Desloc.


n_tot_obj = 4;
global mo  obj_type i_obj
mo = zeros(1,n_tot_obj*2);
obj_type=cell(1,n_tot_obj*2);

main_object=4;% energy
arq.tobj=0;%Multi_Obj
mo([main_object,4+main_object])=1;
 obj_type{main_object}='single';
 obj_type{4+main_object}='single';
 %i_obj = 661; %specific
 arq.tobj=1;%energy
 %arq.tobj=2;%energy
 %arq.tobj=6;%max stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constraint
global ResType n_sigma sign_const

        %arq.tobj=2;%    1---VOLUME  2----ENERGIA   3--Desl.Tot
        %4--Desl.Maxim   5--Desl.Esp   6--Tes.Max
        %arq.tpres=1;%    1---PESO    2----TENSAO     3---Tens e Desloc
        %4---Desloc.

%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy

arq.tpres=[0 0 0 0 0 0];% Reliability

%arq.tpres(1)=1;%volume

arq.tpres(4)=1;%stress
arq.cm = -500;
arq.em =  400;
%  iconst = 3;%stress
%  arq.tpres(iconst)=1;
%  ResType{iconst}='probability';
%  n_sigma(iconst)=3;
%  sign_const(iconst)=1;
%  arq.cm = -5000;
%  arq.em =  5000;
%  
% 
% iconst = 1;%vol
% arq.tpres(iconst)=1;
% ResType{iconst}='deterministic';
% n_sigma(iconst)=0;
% sign_const(iconst)=1;
% % arq.vol0 = 7500;

%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarDist (6,nrv) = [ distribuction; mean; std; par3; par4;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

nRandVar = 14;
var_type(1:6) = 2; % Sensibility variable E
var_type(7:8) = 3; % Sensibility variable Px
var_type(9:14) = 1; % Sensibility variable x



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


RandVarProp(1,:)=var_type; % type var
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
RV_consider = [1:14];
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


%SIGu=arq.MatCurv(end)*.5;