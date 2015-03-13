

%clc,
clear variables global
%clear

global mFext SIGu
global  es0
global  props 
global  offline_stage grad_calculation u_initial

global savedState test_importance
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;




%path('bin',path)
path('HYPLASgui_matlab',path)
path('Mfiles-NonL-POD',path)
path('Stochstcs',path)

grad_calculation=0;
offline_stage=1;
u_initial=[];
test_importance=0;

%path(pathdef);

%% STUDY CASE
%NewPortcLarg2F_RobustProblem
%NewPortcLarg2F_ReliTypeProblem
%NewPortcLarg2F_RRDOTypeProblem
% 
% SpaceTruss_25bar_RobustProblem
%
% Retrieve HYPLAS results file name
%
%HYPLAS_13_6_problem
HYPLAS_kubrick_platewhole

global echo GradOutput_on
echo = 0;
GradOutput_on = 1;
%Paper_TrecLin_RobustProblem

%Eloy10bar_ReliabProb.m

path('Mfmincon',path)

%[Fs,Us,Sigs,ien,vol  ] = Hyplas_analysis( arq.x0, hyplasprojname);
HyplasNLOpt
%TrecNLOpt

rmpath('Mfmincon')

%close; resp;