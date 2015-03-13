

%clc,
clear variables global
%clear

global mFext var_type ihistogram SIGu test_importance
global  es0 sens NLGEOM HYPLAS_flag
global  props echo echoReliab echoReliabPlot GradOutput_on
global  offline_stage grad_calculation u_initial

%% Setting Options
GradOutput_on = 1;
echo = 0;
echoReliab = 0;
echoReliabPlot=0;
sens = 1;
grad_calculation=0;% NL_fesol 
u_initial=[];

offline_stage=0; % Ensure new Z
ioffT = 5;
noffT = 60;

test_importance=0;
ihistogram = 0;
var_type = 1; % Sensibility variable x


NLGEOM = 0;
disp('NLGEOM')
disp([NLGEOM])

%FEM ANALYSIS
%arq=arqfileGeneret();
%load ExempA-Messac.mat
%load 200barLeoRB.mat
% load arq0
% load arq4bar
% arq.nox
% global Messac
% Messac=0;

%% Setting Stream seed (random generator)
global savedState 
Rebuilt_Stream = 0;
if Rebuilt_Stream 
    defaultStream = RandStream.getDefaultStream;
    save Arq_saved_stream defaultStream
else
    load Arq_saved_stream
end
savedState = defaultStream.State;

%% Add Paths
path(pathdef);
path([pwd '\HYPLASmat_files'],path)
path([pwd '\Examples'],path)
path([pwd '\Mfmincon'],path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STUDY CASE
% NewPortcLarg2F_RobustProblem
%   NewPortcLarg2F_ReliTypeProblem

% PAPER RMO_POD_RBRDO %
%   NewPortcLarg2F_RRDOTypeProblem
%   NewPortcLarg2F_RRDOTypeProblem2
%    NewPortcLarg2F_RRDOTypeProblem2obj

% SpaceTruss_25bar_RobustProblem
% SpaceTruss_25bar_RobustProblemR2BDO

%SpaceTruss_Dome_RobustProblem

%Paper_TrecLin_RobustProblem

%Eloy10bar_ReliabProb.m
% path([pwd '\Mfiles-NonL-POD'],path)

HYPLAS_kubrick_platewhole
% KubrickSqrtHole_PerfcPlas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RBtype
    ioffT = N; % No convergence test
    noffT = N; % No convergence test
end
global  mat_curv linear_prob
% COLAR NOS ARQUIVOS DE ENTRADA
%arq=MatProp(arq);
%SIGu=mean(arq.MatCurv(3:4,2));

if HYPLAS_flag
    cd('..')
    bindir = [pwd '\' bindir];
    cd('Mfiles-NonL-POD')
    HYPLASexe='hyplas90a_OK.exe';
    funame = 'Hyplas_analysis';
else
    mat_curv=arq.MatCurv;
    % LINEAR PROBLEM
    if linear_prob
        esc=1e50;
        arq.MatCurv=[esc/arq.E(1), esc];
        mat_curv=arq.MatCurv;
    end

    [arq,area,props,fext,glb,tipo,nel,mdf,mdv]=pre_compus_arq(arq);
    fext =fext*k_p;

    mFext=max(abs(fext));

    es0=zeros(nel,1);
end

%%%Simple FEM Analysis
%[Fs,Us,Sigs,arq.ien,arq.ivol]=NLfesol(area,props,fext,glb,tipo);

TrecNLOpt

rmpath([pwd '\MO'])

%close; resp;