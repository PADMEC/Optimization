function varargout=post_proc(Resum,iplot)
if nargin<2
    iplot=1;
end
% Objectives
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd
Objectives={'Volume'; 'Displ'; '\sigma'; 'Energy'; '\lambda'};
iftc=0;ft=[];
nrtot=length(Resum);
for i=1:nrtot
    % Read Results (Resum)
    if length(Resum{i})==7;
        [f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp]=Resum{i}{:};
        
    elseif length(Resum{i})==8;
        [f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp,mo]=Resum{i}{:};
    else
        fprintf('Results %d do not match, %d \n',i,length(Resum{i}))
        continue
    end
    
    % Preper global solution
    %ifft=[1 36*(1:4)+1;36*(1:5)];
    iftc=iftc+1;
    ifft(1,i)=size(ft,1)+1;
    ft=[ft;f];
    ifft(2,i)=size(ft,1);
end

% Find geral domained solutions (ft)
[pnondom,whodom]=filtr_2_pdom(ft);
if isempty(whodom), pdom=whodom;
else pdom=whodom(:,1)';
end

fprintf('Results MO methods: Time, F Count, Num. NonPareto, Evness, Pareto Area \n')
for i=1:nrtot
    Objectivesi=Objectives;
    % Read Results (Resum)
    if length(Resum{i})==7;
        [f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp]=Resum{i}{:};
        objs=[];
        
    elseif length(Resum{i})==8;
        [f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp,mo]=Resum{i}{:};
        objs=find(mo);
        nmo=length(mo);
        if nmo/2>=4
            % Statistics objectives (RO)
            for iobj=1:nmo/2
                Objectivesi{iobj+nmo/2}=['SD ' Objectives{iobj}];
                Objectivesi{iobj}=['Mean ' Objectives{iobj}];
            end
        end
    else
        continue
    end
    
    Tot_Fcount=sum(fcount);

    % Non Converged Solutions
    nonConv=find(Cnvrg<0);
    
    % Separate domined for solutions i
    pdomtot=pdom((pdom>=ifft(1,i))&(pdom<=ifft(2,i)));
    pdomi=pdomtot-ifft(1,i)+1;
%     [pdom]=filtr_2_pdom(f,f);
    nonPP=union(nonConv,pdomi);%[];%
    nnpi=length(nonPP);
    fok=f;
    fok(nonPP,:)=[];
    ferr=f(nonPP,:);
    
    %PLOTS
    if iplot
        %ierotim=find(Cnvrg<1);
%         %iok=find(Cnvrg>0);
%         fok=f(Cnvrg>=0,:);
%         ferr=f(Cnvrg<0,:);
%         if size(f,2)>2
            Pareto_Plot(fok,ferr,{Objectivesi{objs}},MOptype)
%         end
    end

    ndp=length(Cnvrg);
    %  MO_Geom_ProbGraf

    % Evenness (ev) parameter computation
    [ev,dl,du,Aregion,tr]=UnifmComput(fok,ptsp);

    iok=1:ndp;
    %m0=find(du<10^-7);
    iok(nonPP)=[];iok(tr)=[];
    %iok(du<10^-6)=[];
    ndpP(i)=length(iok);
    %post_proc(f(iok,:),x,Cnvrg(iok),fcount,Tot_Time,MOptype,ptsp);

    TTs(i)=Tot_Time;
    TFCs(i)=Tot_Fcount;
    eveness(i)=ev;
    Par_area1(i)=sum([dl;du])/2/ev;
    Par_area(i)=Aregion;
    nnonP(i)=nnpi;
    fprintf('%d, %d, %d, %d, %d \n',Tot_Time, Tot_Fcount, nnpi, ev, Par_area(i))
end
if iplot
    ftok=ft(pnondom,:);
    fter=ft(pdom,:);
    Pareto_Plot(ftok,fter,{Objectivesi{objs}},'MO Filtered Pareto Solutions')
end
varargout={TTs,TFCs,eveness,Par_area1,Par_area,nnonP,ndpP};
