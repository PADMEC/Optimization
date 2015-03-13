%SpaceTruss_Dome

arq.name = 'SpaceTruss_Dome';

nodes=[      0         0    8.2160
       25.0000         0    6.2160
       12.5000   21.6500    6.2160
      -12.5000   21.6500    6.2160
      -25.0000         0    6.2160
      -12.5000  -21.6500    6.2160
       12.5000  -21.6500    6.2160
       43.3000  -25.0000         0
       43.3000   25.0000         0
             0   50.0000         0
      -43.3000   25.0000         0
      -43.3000  -25.0000         0
             0  -50.0000         0];

% c1=[ones(6,1),(2:7)'];
% c2=[(2:7)',[3:7,2]'];
% g=[8:13,8];p=[2:7,2];
% aa(1:2:14)=g;
% aa(2:2:14)=p;
% c3=[[7 aa]',[aa 8]'];
% connect =[c1;c2;c3];
connect =[1     2
     1     3
     1     4
     1     5
     1     6
     1     7
     2     3
     3     4
     4     5
     5     6
     6     7
     7     2
     7     8
     8     2
     2     9
     9     3
     3    10
    10     4
     4    11
    11     5
     5    12
    12     6
     6    13
    13     7
     7     8
     8     2
     2     8];

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
    grid on
end

% Loads
%initializing
arq.px=zeros(1,length(arq.nox));
arq.py=zeros(1,length(arq.nox));
arq.pz=zeros(1,length(arq.nox));
% setting
arq.pz(1)=-20; %P1
%arq.pz([2:7])=-1; %P2
arq.sdpz(1)=1;

%Fix nodes
arq.ax=8:13;
arq.ay=8:13;
arq.az=8:13;

%Material Properties
arq.ro=0.1;

arq.E=zeros(1,arq.nb);
arq.E(1:end)=21000;

%Links
arq.arex=ones(arq.nb,1);
arq.arex( 1: 6,1)=1;
arq.arex( 7:12,1)=2;
arq.arex(13:27,1)=3;
% nvs=[6,6,15];
% currv=1;
% for iv=1:length(nvs)
%     endv=currv+nvs(iv)-1;
%     disp([currv,endv])
%     arq.arex(currv:endv,1)=iv;
%     currv=currv+nvs(iv);
% end

%Optimization Settings
arq.x0=ones(1,6);
arq.xmin=0.1;
arq.xmax=1000;

arq.vol0=750;
arq.cm=-5000;
arq.em=5000;


%Random Variables
arq.sdE=zeros(length(connect));
arq.sdpx=zeros(1,length(arq.nox));
arq.sdpy=zeros(1,length(arq.nox));
arq.sdpz=zeros(1,length(arq.nox));
arq.sdx=zeros(1,6);

%arq.sdE(1)=2e5./arq.E(1);
