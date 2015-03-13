
global G_const G_limits betatarg RBtype linStressComputed linear_prob
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl
global serieLimitState Stochast_moment GradOutput_on DinamcVar Method
linStressComputed=0;
GradOutput_on=1;

Stochast_on=1;

SpaceTruss_Dome
 
arq.x0=[1 1 1]*10;
 
%x([2:end,1],:),%NBI
%    45.5439    2.1181    0.5163
%    45.7016    3.4252    0.7578
%    44.0692    4.2809    1.0052
%    42.1487    5.4306    1.2511
%    40.4164    6.3474    1.5109
%    37.9828    8.0555    1.7430
%    36.3809    8.7988    2.0163
%    33.2417   11.2442    2.2390
%    30.6843   13.0071    2.4930
%    27.6319   15.2598    2.7492
    
%arq.x0=[1 1 1]*20;
  X_ops=arq.x0;


%RV_consider = [6,11,13,14];
%DinamcVar = 3;Method='PC';


 
%N/mm^2 = MPa = 1MN/(1e2cm)^2 = 1MN/1e4cm^2 = 1e2N/cm^2
 %2e5MPa = 2e7 N/cm^2

Stochast_moment = 1;
%DinamcVar = 10000;Method='MC';
DinamcVar = 2;Method='PC';
%DinamcVar = 2;Method='R2DO';

%DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';


%k_p = 2000;%500
k_p = 1;
%RBtype=1;%POD
%RBtype=2;%RBM
RBtype=0;%FEM
linear_prob=1;
%RBtype=0;
 
 %% Constraint
global ResType n_sigma sign_const specfc_FF MinLambd

%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd

arq.tpres=[0 0 0 0 0 0];% Reliability

 iconst = 2;%displ
 arq.tpres(iconst)=1;
 ResType{iconst}='probability';
 n_sigma(iconst)=3;
 sign_const(iconst)=1;
 arq.d = [0 0 1.];
 specfc_FF(iconst) = 3; %
  
 iconst = 5;%critc load
 arq.tpres(iconst)=1;
 ResType{iconst}='probability';
 n_sigma(iconst)=-5;
 sign_const(iconst)=1;%Maximize
 specfc_FF(iconst) = 0; %
 MinLambd = 1;

iconst = 1;%vol
arq.tpres(iconst)=1;
ResType{iconst}='deterministic';
n_sigma(iconst)=0;
sign_const(iconst)=1;
arq.vol0 = 775;

%% Objective
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5-Lambd

n_tot_obj = 5;
global mo  obj_type sign_obj
mo = zeros(1,n_tot_obj*2);
obj_type=cell(1,n_tot_obj*2);

%main_object=3;% energy
main_object=5;% 
arq.tobj=0;%Multi_Obj
mo([main_object,n_tot_obj+main_object])=1;
 obj_type{main_object}='single';
 obj_type{main_object+n_tot_obj}='single';
 
 %Vol
%  mo(1)=1;
%  obj_type{1}='single';
 
 % minimization
 sign_obj = ones(1,sum(mo));
 if mo(5)
     % maximization of critical load
     sign_obj(sum(mo(1:5)))=-1;
 end
 %i_obj = 661; %specific
 %arq.tobj=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarDist (6,nrv) = [ distribuction; mean; std; par3; par4;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

nRandVar = 5;
randgrp={1:3,4,5};
var_tp(randgrp{1}) = 1; % Sensibility variable x
var_tp(randgrp{2}) = 2; % Sensibility variable Ex
var_tp(randgrp{3}) = 3; % Sensibility variable P

E=arq.E;dE=E/5;% 5,0%


%Random Variables
% 1 - E all design regions (3)
% 2 - Load (1)
RandVarProp(1,:)=var_tp; % type var
RandVarProp(2,randgrp{1})=1; % which type data
RandVarProp(2,randgrp{2})=1; % which type data 
RandVarProp(2,randgrp{3})=3; % which type data

%(DoF node i direction x = i*3-2)
% Components
rv_data_comp = {1,2,3,1:3,3};

RandVarDist(1,:)=2*ones(1,nRandVar); % distribuction: 2-LN
RandVarDist(2,randgrp{1})=20; % Mu par1
RandVarDist(2,randgrp{2})=21000; % Mu par1
RandVarDist(2,randgrp{3})=20; % Mu par1

RandVarDist(3,randgrp{1})=1; % std par2
RandVarDist(3,randgrp{2})=1050; % std par2
RandVarDist(3,randgrp{3})=3; % std par2
%RandVarDist(3,9:14)=0.05; % std par2

%%%%%%%%%%%%%%%%%%%%%%%%
% R. V. Elimination!! 
% R. V. Considered:
RV_consider = 1:nRandVar;
% Re-setting random var data
nRandVar = length(RV_consider);
rv_data_comp = {rv_data_comp{RV_consider}};
RandVarDist=RandVarDist(:,RV_consider);
RandVarProp=RandVarProp(:,RV_consider);
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit State (G=0) CONSTRAINT
G_fun = 'StructLimitState';
S_fun = 'Stoch_Analysis';
%G_const = 1 VOLUME, 2 FLAMBAGEM GLOBAL, 3 FLAMBAGEM LOCAL, 4 TENSAO, 
% 5 DESLOCAMENTO, 7 DESLOCAMENTO ESPECIFICO, 8 TENSAO ESPECIFICA

% NEW!!!
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5-Lambd


%G_const = [3,4];
%G_limits = [1,10e3];
% G_const = 4;
% G_limits = 65000;

 G_const = 5;
 G_limits = -10;
 
serieLimitState = [];

betatarg = 3.;
%Mu = [1,1,1]*.5;
%Cu = diag([1 2 3]/50);
%ft=[2,2,2];


nrv=size(RandVarProp,2);
Correl=eye(nrv);

arq.modos = [1,1,0];
%arq=MatProp(arq);
%SIGu=arq.MatCurv(end)*.5;