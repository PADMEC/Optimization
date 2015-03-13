
global G_const G_limits betatarg specfc_FF ... 
RBtype linStressComputed linear_prob

global Stochast_on RandVarProp rv_data_comp ...
    RandVarDist Correl DinamcVar Method

global serieLimitState Stochast_moment
linStressComputed=0;
Stochast_on=1;

% Created by Renato S. Motta (03/2014)
global bindir HYPLASexe
global NO_IG HYPLAS_flag hyplasprojname
NO_IG=1;
HYPLAS_flag=1;
bindir = 'bin';
HYPLASexe = 'dhyplas90a.exe';

hyplasprojname='KubrickSqrtHole';

% x = [250,250];

% 
 x = [120,100];
% 
 X_ops = [
%   357.6893  276.1824
%   391.3865  100.0000
   100.0000  100.0000];
% x=X_ops(1,:);

arq.x0=x;

arq.xmin=100;
arq.xmax=477;

nvar = length(arq.x0);
arq.arex = 1:nvar;

Stochast_moment = 1;
%  DinamcVar = 10000,Method='MC'
 DinamcVar = 2,Method='PC'
%  DinamcVar = 2;Method='R2DO';

%DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';

%  Method='R2DO' , %PC method with FORM
%  Method='R2DOmc' , %PC method with MC

%k_p = 2000;%500
k_p = 1;
 RBtype=1;N=30;%POD1
%RBtype=2;%RBM2
% RBtype=0;%FEM0
linear_prob=0;
%RBtype=0;

%% Objective
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};

        %arq.tobj=2;%    1---VOLUME  2----ENERGIA   3--Desl.Tot
        %4--Desl.Maxim   5--Desl.Esp   6--Tes.Max
        %arq.tpres=1;%    1---PESO    2----TENSAO     3---Tens e Desloc
        %4---Desloc.

%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy

n_tot_obj = 4;
global mo  obj_type i_obj
mo = zeros(1,n_tot_obj*2);
obj_type=cell(1,n_tot_obj*2);

arq.tobj=0;%0-Multi_Obj, n-single objective
% main_object=4;% energy
main_object=3;% stress
mo([main_object,n_tot_obj+main_object])=1;
  obj_type{main_object}='specific';
  obj_type{4+main_object}='specific';
  i_obj([1:2]) = [1,1]*264; %specific [0,120;8.5605,119.6930]

% main_object=2;% displacement
% mo([main_object])=1;
% obj_type{main_object}='specific';
% %   i_obj(3) = 3; %displacement node 2_x
%   i_obj(3) = 5; %displacement node 3_x

% % Third objective (MO3)
% main_object=1;% displacement
% mo([main_object])=1;
% obj_type{main_object}='max';
  
%  arq.tobj=1;%vol
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

Volum_constraint=0;
if Volum_constraint
% arq.tpres(1)=1;%volume
     iconst = 1;%volume
     arq.tpres(iconst)=1;
     ResType{iconst}='deterministic';
     n_sigma(iconst)=0;
     sign_const(iconst)=1;
end
if  ~strcmp(Method,'R2DO')
    % No reliabilit analysis
     iconst = 2;%displ
     arq.tpres(iconst)=1;
     ResType{iconst}='probability';
     n_sigma(iconst)=3.2;
     sign_const(iconst)=1;
     arq.d(8) = 6;
     specfc_FF(iconst) = 8;% y displacement of the 4th node
%      specfc_FF(iconst) = 3;% x displacement of the 2th node
end
 
%  
% 
% iconst = 1;%vol
% arq.tpres(iconst)=1;
% ResType{iconst}='deterministic';
% n_sigma(iconst)=0;
% sign_const(iconst)=1;
% arq.vol0 = 7500;

%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarDist (6,nrv) = [ distribuction; mean; std; par3; par4;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

%nRandVar = 14;
% var_type(1:6) = 2; % Sensibility variable E
% var_type(7:8) = 3; % Sensibility variable Px
% var_type(9:14) = 1; % Sensibility variable x

%Random Variables
% 1 - E all design regions (3)
% 2 - Load (1)
% RandVarProp(1,:)=[2 2 2 3]; % type var
% RandVarProp(2,:)=[1 1 1 5]; % which type data
% rv_data_comp = {1,2,3,0};
% 

RandVarDist(1,1:3)=[1,1,2]; % distribuction
RandVarDist(2,1:3)=[0,0,1]; % Mu par1 (F, +40, +80)
RandVarDist(3,1:3)=[.5,.2,.2]; % std par2 (F, *40 = 20,8)

%%%%%%%%%%%%%%%%%%%%%%%%
%%R. V. Elimination!! 
%%R. V. Considered:
% nRandVar = size(RandVarDist,2);
% RV_consider = [1:nRandVar];
 RV_consider = [1:2];
nRandVar = length(RV_consider);
%rv_data_comp = {rv_data_comp{RV_consider}};
RandVarDist=RandVarDist(:,RV_consider);
%RandVarProp=RandVarProp(:,RV_consider);
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit State (G=0)
% G_fun = 'StructLimitState';
% S_fun = 'Stoch_Analysis';
G_fun = 'StructLimitState';
% G_fun = 'HyplasLimitState';
S_fun = 'Stoch_HyplasAnalysis';

%G_const = 1 VOLUME, 2 DESLOCAMENTO, 3 TENSAO, 4 COMPILANCE
%vol,u,vonmis,sn
%G_const = [3,4];
%G_limits = [1,10e3];
% G_const = 4;
% G_limits = 65000;

global set_Failure_FORM
if strcmp(Method,'R2DO')
    G_const = 2;% displacement
    G_limits = 6;% maximum displacement value
    set_Failure_FORM = 8;% y displacement of the 4th node
    serieLimitState = [];
    betatarg = 3.2;
end
%Mu = [1,1,1]*.5;
%Cu = diag([1 2 3]/50);
%ft=[2,2,2];


nrv=size(RandVarDist,2);
Correl=eye(nrv);


%SIGu=arq.MatCurv(end)*.5;