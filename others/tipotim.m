function arq=tipotim(arq)


obj=findobj('Tag','fpes');
if length(obj)>1
    parnt=get(obj(2:end),'parent');
    parnt=get(parnt,'parent');
    close(parnt)
    obj(2:end)=[];
end
fpes=get(obj,'value');

obj=findobj('Tag','fen');
fen=get(obj,'value');

obj=findobj('Tag','fflamb');
flamb=get(obj,'value');

obj=findobj('Tag','fdes');
fdes=get(obj,'value');

obj=findobj('Tag','fdesesp');
fdesesp=get(obj,'value');

obj=findobj('Tag','fdest');
fdest=get(obj,'value');

obj=findobj('Tag','ftes');
ftes=get(obj,'value');

fde=0;fte=0;
if fdesesp==1
    if isfield(arq,'idm')
        fde=length(arq.idm);
    end
    if isfield(arq,'kb')
        fte=length(arq.kb);
    end
    fdesesp=fde+fte;
else
    if isfield(arq,'idm')
        arq=rmfield(arq,'idm');
    end
    if isfield(arq,'kb')
        arq=rmfield(arq,'kb');
    end
end
sobj=fpes+fen+flamb+fdest+fdes+ftes+fdesesp;
arq.modos=[1 flamb 0];%vib
if sobj==1;
    switch 1
    case fpes
        arq.tobj=1;
    case fen
        arq.tobj=2;
    case flamb
        arq.tobj=3;
    case fdest
        arq.tobj=4;
    case fdes
        arq.tobj=5;
    case ftes
        arq.tobj=6;        
    case fdesesp
        arq.tobj=7;
    otherwise
        errrr_fobj_tipotim
    end
else
    arq.tobj=0;
    arq.mo=[fpes fen flamb fdest fdes ftes fdesesp];
end
obj=findobj('Tag','resdes');
frd=get(obj,'value');

obj=findobj('Tag','resinsglb');
frig=get(obj,'value');

obj=findobj('Tag','resinsloc');
fril=get(obj,'value');

obj=findobj('Tag','restes');
frt=get(obj,'value');

obj=findobj('Tag','respes');
frp=get(obj,'value');

arq.tpres=[frp frig fril frt frd];

%%%%%%%%%Restrição de flambagem global, modo de flamb = 1
if frig>0
    arq.modos(2)=1;
end

%%%%%%%%%Restrição de frequencia, modo de vibr = 1
%if frfreq>0
%    arq.modos(3)=1;
%end

%arq.tpres=[frp frt frd];