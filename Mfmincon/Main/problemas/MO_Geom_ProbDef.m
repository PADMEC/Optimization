%Geometric problem definition.

%Objectives
optPoints=[1,0;
           0,1
          -1,0]%
          %0,-1];
%Constraints
Rc=.5;
Cc=[1 .4];

global funinput coninput
nOptPt=3;nRestC=1;
funinput={optPoints,nOptPt};
coninput={Cc,Rc,nRestC};