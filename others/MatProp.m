function arq=MatProp(arq)


%ANALYSIS OPTIONS
%arq.modos(1) - Elastc Analysis
%arq.modos(2) - Buckling Analysis
%arq.modos(3) - Vibration Analysis
arq.modos= [1 0 0];
nelm=size(arq.barra,1);
%Verificar UNIDADES
    is=-16:16;OG=10.^(is);

    m=max([arq.nox(:); arq.noy(:)]);
    iOG=find(m<OG,1)-1;
    m_dim=OG(iOG);

    f=max(abs([arq.px arq.py]));
    iOG=find(f<OG,1)-1;
    f_dim=OG(iOG);
    
%arq.name= 'Flamb3X.mat';
if ~isfield(arq,'E')

    E= 2.0000e+08;%KN/m2=2e5MPa
    arq.E=E*m_dim^2/f_dim;
end
if ~isfield(arq,'x0')
    x0= 5e-4/nelm;%m2
    arq.x0=x0*m_dim^2;
end


%arq.x0= 1.0000e-003;
%arq.E= 0;
E=arq.E(1);
esc=500*1e3;%KN/m2=500MPa
esc=esc*E/2.0000e+08;
arq.MatCurv=[esc/E, esc
            esc/E*1.5, esc*1.2
            esc/E*2, esc*1.4
            esc/E*4, esc*1.6
            esc/E*7, esc*1.8
            esc/E*14, esc*2.
            esc/E*30, esc*2.2
            esc/E*200, esc*2.5];


%arq.xmin= 1.0000e-006*m_dim^2;
%arq.xmax= 0.0100*m_dim^2;
%arq.ro= 7.85*m_dim^3/f_dim;
arq.nfh= 20;
arq.nfv= 20;

%Obje Function
%1-pes; 2-en; 3-flamb; 4-dest; 5-des; 6-tens; 7-des esp;
%arq.tobj= 1;
%Constrn
%arq.tpres= [0 0 1 1 1];
arq.Converg= 1;
arq.mo= [0 0 0 0 0 0];%mult obj optim - obj functions
%arq.cm= -500000000;
%arq.em= 500000000;
%arq.d= 0.0200*m_dim;
