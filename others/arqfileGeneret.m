function arq=arqfileGeneret()

coords=[
  -0.5000    0.5000
         0    0.5000
   -0.5000    1.0000
         0    1.0000
   -0.5000    1.5000
         0    1.5000
   -0.5000    2.0000
         0    2.0000
   -0.5000    2.5000
         0    2.5000
   -0.5000    3.0000
         0    3.0000
   -0.5000    3.5000
         0    3.5000
   -0.5000    4.0000
         0    4.0000
   -0.5000    4.5000
         0    4.5000
   -0.2500         0
   -0.2500    5.0000];


conect=[
     1     3
     1     4
     1     2
     2     3
     2     4
     3     5
     3     6
     3     4
     4     5
     4     6
     5     7
     5     8
     5     6
     6     7
     6     8
     7     9
     7    10
     7     8
     8     9
     8    10
     9    11
     9    12
     9    10
    10    11
    10    12
    11    13
    11    14
    11    12
    12    13
    12    14
    13    15
    13    16
    13    14
    14    15
    14    16
    15    17
    15    18
    15    16
    16    17
    16    18
    17    18
    19     1
    19     2
    17    20
    20    18];


Link=[  2, 3, 2, 3, 2, 2, 3, 3, 3, 2, 2, 3, 3, 3, 2, 2, 3, 3, 3, 2, 2,...
  3, 3, 3, 2, 2, 3, 3, 3, 2, 2, 3, 3, 3, 2, 2, 3, 3, 3, 2, 2, 1, 1, 1, 1];

%Coords
arq.nox= coords(:,1)';
arq.noy= coords(:,2)';
%LOADs
arq.px= [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
arq.py= [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.3];
%Conections
arq.barra= conect;
%Bond Cond Dircl
arq.ax= [19 20]; 
arq.ay= 19;

%ANALYSIS OPTIONS
%arq.modos(1) - Elastc Analysis
%arq.modos(2) - Buckling Analysis
%arq.modos(3) - Vibration Analysis
arq.modos= [1 0 0];

arq.name= 'Flamb3X.mat';
arq.x0= 1.0000e-003;
arq.E= 2.0000e+05;
%arq.E= 0;
E=2.0000e+05;
esc=700;
arq.MatCurv=[esc/E, esc
            esc/E*2, esc*1.1]


global  mat_curv
mat_curv=arq.MatCurv;
        
arq.xmin= 1.0000e-006;
arq.xmax= 0.0100;
arq.ro= 7850;
arq.arex= Link;
arq.nfh= 20;
arq.nfv= 20;
arq.nb= 45;

%Obje Function
%1-pes; 2-en; 3-flamb; 4-dest; 5-des; 6-tens; 7-des esp;
arq.tobj= 1;
%Constrn
arq.tpres= [0 0 1 1 1];
arq.Converg= 1;
arq.mo= [0 0 0 0 0 0];%mult obj optim - obj functions
arq.cm= -500000000;
arq.em= 500000000;
arq.d= 0.0200;
