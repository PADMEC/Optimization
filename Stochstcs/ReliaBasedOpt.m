

clear variables global

global savedState echo
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;

grad_calculation=0;
u_initial=[];
echo=0;


global G_const G_limits betatarg RBtype linStressComputed
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl DinamcVar Method
global serieLimitState Stochast_moment var_type GradOutput_on
linStressComputed=0;
GradOutput_on=1;

 
 %% Constraint
global ResType n_sigma sign_const

% iconst = 3;%disp
% arq.tpres(iconst)=1;
% ResType{iconst}='probability';
% n_sigma(iconst)=3;
% sign_const(iconst)=1;
% arq.d = 2;

iconst = 1;
arq.tpres(iconst)=1;
ResType{iconst}='deterministic';
n_sigma(iconst)=0;
sign_const(iconst)=1;
%% Objective

f = 'randF';
Stochast_on=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor
%RandVarProp(1,:)=[3 3]; % type var
%RandVarProp(2,:)=[4 4]; % which type data
%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarDist (6,nrv) = [ distribuction; par1; par2;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

RandVarDist(1,:)=[1 1 ]; % distribuction
RandVarDist(2,:)=[4 4]; % Mu par1
RandVarDist(3,:)=[1 1]; % std par2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit State (G=0)
%G_fun = 'StructLimitState';
%S_fun = 'Stoch_Analysis';
%G_const = 1 VOLUME, 2 FLAMBAGEM GLOBAL, 3 FLAMBAGEM LOCAL, 4 TENSAO, 
% 5 DESLOCAMENTO, 7 DESLOCAMENTO ESPECIFICO, 8 TENSAO ESPECIFICA

 G_const = 5;
 G_limits = 1.2;
 
serieLimitState = [];

betatarg = 3.;
%Mu = [1,1,1]*.5;
%Cu = diag([1 2 3]/50);
%ft=[2,2,2];

nrv=size(RandVarDist,2);
Correl=eye(nrv);

t_ord = 0;

x0=25;

%Stochast_moment = 1;
%DinamcVar = 1000;Method='MC';
Stochast_moment = 0;
DinamcVar = 0;Method='FORM';
i=0;
%for x0=.5:.1:2
    [PF,beta, dbeta]=FunStatistcalc(f, x0,t_ord);
    i=i+1;Res(i,:)=[PF,beta, dbeta];
%end
%return
Method='MC';
n_Conv_MC=0;
n_x_var=20;
DinamcVar = ceil(10/PF)
DinamcVar = 1000;
[PF,b,db]=FunStatistcalc(f, x0,t_ord);

%DinamcVar = 0;Method='FORM';
if n_x_var
    xs=x0-n_x_var/2 + [1:.2:n_x_var];
    for im=1:length(xs)
        xi = xs(im);
        [PF,b,db]=FunStatistcalc(f, xi,t_ord);
        disp([PF*100,im])
        pfs(im)=PF;
        bs(im)=b;
    end
    figure,
    plot(xs,pfs,'.-b')
end

if n_Conv_MC
    Ns=[5:5:50].^3;
    for im=1:30
        for in=1:length(Ns)
            DinamcVar = Ns(in);
            [PF,b,db]=FunStatistcalc(f, x0,t_ord);
            disp([PF*100,in])
            pfs(im,in)=PF;
            bs(im,in)=b;
        end
    end
    figure,
    plot(Ns'*ones(1,30),pfs','.-b')
    hold on, plot(Ns,mean(pfs),'--k')
end

return
Stochast_moment = 1;
DinamcVar = 5;Method='PC';
[ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ]=FunStatistcalc(f, x0);
-M_out{1}/S_out{1}