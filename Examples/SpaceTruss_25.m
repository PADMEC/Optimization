%SpaceTruss_25

arq.name = 'SpaceTruss_25';

nodes = [  -37.5000         0  200.0000
   37.5000         0  200.0000
  -37.5000   37.5000  100.0000
   37.5000   37.5000  100.0000
   37.5000  -37.5000  100.0000
  -37.5000  -37.5000  100.0000
 -100.0000  100.0000         0
  100.0000  100.0000         0
  100.0000 -100.0000         0
 -100.0000 -100.0000         0];

connect =[     1     2
     2     3
     1     4
     2     6
     1     5
     2     4
     2     5
     1     3
     1     6
     3     6
     5     4
     3     4
     5     6
    10     3
     6     7
     4     9
     5     8
     7     4
     3     8
    10     5
     6     9
    10     6
     7     3
     4     8
     5     9];

arq.nox=nodes(:,1);
arq.noy=nodes(:,2);
arq.noz=nodes(:,3);

arq.barra=connect;
arq.nb=length(arq.barra);

lplot=0;
if lplot
    %xd=arq.nox;xd(glx)=xd();
    figure,plot3([arq.nox(connect(:,1)),arq.nox(connect(:,2))]',...
    [arq.noy(connect(:,1)),arq.noy(connect(:,2))]',...
    [arq.noz(connect(:,1)),arq.noz(connect(:,2))]','.-')
end

% Loads
arq.px=zeros(1,length(arq.nox));
arq.py=zeros(1,length(arq.nox));
arq.pz=zeros(1,length(arq.nox));
arq.pz([1,2])=-1e4;
arq.py([1,2])=-1e4;
arq.px([3,6])=500;
arq.sdpx([3,6])=50;

%Fix nodes
arq.ax=7:10;
arq.ay=7:10;
arq.az=7:10;

%Material Properties
arq.ro=0.1;

arq.E=zeros(1,arq.nb);
arq.E(1)=1e7;
arq.E(2:5)=1e7;
arq.E(6:9)=1e7;
arq.E(10:13)=1e7;
arq.E(14:21)=1e7;
arq.E(22:25)=1e7;

%Links
arq.arex=ones(length(connect),1);
arq.arex(1,1)=1;
arq.arex(2:5,1)=2;
arq.arex(6:9,1)=3;
arq.arex(10:13,1)=4;
arq.arex(14:21,1)=5;
arq.arex(22:25,1)=6;

%Optimization Settings
arq.x0=ones(1,6);
arq.xmin=0.05;
arq.xmax=10;

arq.vol0=750;
arq.cm=-5000;
arq.em=5000;


%Random Variables
arq.sdE=zeros(length(connect));
arq.sdpx=zeros(1,length(arq.nox));
arq.sdpy=zeros(1,length(arq.nox));
arq.sdpz=zeros(1,length(arq.nox));
arq.sdx=zeros(1,6);

arq.sdE(1)=2e5./arq.E(1);
arq.sdE(2:5)=2e5./arq.E(2:5);
arq.sdE(6:9)=2e5./arq.E(6:9);
arq.sdE(10:13)=2e5./arq.E(10:13);
arq.sdE(14:21)=2e5./arq.E(14:21);
arq.sdE(22:25)=1.5e6./arq.E(22:25);

arq.sdpx(3)=50;
arq.sdpx(6)=50;

arq.sdx(1:6)=0.05;
