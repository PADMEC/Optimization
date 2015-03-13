
global G_const G_limits betatarg RBtype 
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl DinamcVar Method

 %load  NewPortcLarg
 load  10barvar
 %arq.x0=[0.1841    1.0218    1.3823];%Lin FEM
 arq.x0=[ 0.1820    1.0136    1.3914];%Lin POD
 arq.x0=[ 0.2104    1.0396    1.3516];%Lin POD
 %arq.x0=[1 1 1];
 arq.xmin=0.2;
 arq.E=2.07e7;
 %k_p =10000;
 %k_p =6000;
 k_p = 5000;%500
RBtype=1;
arq.tpres=[0 0 0 0 0 1];% Reliability
%arq.tobj=6;%sig
arq.tobj=1;%vol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarProp (6,nrv) = [ distribuction; par1; par2;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

Stochast_on=1;
DinamcVar = 1000;Method='MC';
% DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';
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


RandVarProp(1,:)=[2  3]; % type var
RandVarProp(2,:)=[1  5]; % which type data
rv_data_comp = {[1,2,3],0};

RandVarDist(1,:)=[2  1]; % distribuction
RandVarDist(2,:)=[E  1]; % Mu par1
RandVarDist(3,:)=[dE  .25]; % std par2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit State (G=0)
G_fun = 'StructLimitState';
%G_const = 1 VOLUME, 2 FLAMBAGEM GLOBAL, 3 FLAMBAGEM LOCAL, 4 TENSAO, 
% 5 DESLOCAMENTO, 7 DESLOCAMENTO ESPECIFICO, 8 TENSAO ESPECIFICA
%G_const = [3,4];
%G_limits = [1,10e3];
G_const = 4;
G_limits = 51750;
betatarg = 3.3;
%Mu = [1,1,1]*.5;
%Cu = diag([1 2 3]/50);
%ft=[2,2,2];


nrv=size(RandVarProp,2);
Correl=eye(nrv);