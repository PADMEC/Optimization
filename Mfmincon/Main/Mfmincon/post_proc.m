function varargout=post_proc(f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp)


%PLOTS
if 1==1
%f=[f(:,3) f(:,1) f(:,2) ];

%ierotim=find(Cnvrg<1);
%iok=find(Cnvrg>0);
fok=f(Cnvrg>0,:);
ferr=f(Cnvrg<1,:);


mxf=(max(f));
mnf=(min(f));
%   mxf =[395   395     0]
%   mnf =[ 0     0   -14]
if size(f,2)==3
    figure1=figure;
    %    hold on
        %chin=fill3(ft(:,1),ft(:,2),ft(:,3),'b');
        %alpha(chin,.3);
        %plot3(bt(:,1),bt(:,2),bt(:,3),'.')
        

            axes1 = axes('Parent',figure1,'FontWeight','bold','FontSize',14);
            grid('on');
            hold('all');

            % Create scatter3
            %scatter3(X1,Y1,Z1,S1,C1,'MarkerFaceColor','flat','MarkerEdgeColor','none',...
            %    'Parent',axes1);

        
        scatter3(fok(:,1),fok(:,2),fok(:,3),50,fok(:,3),'filled')
        plot3(ferr(:,1),ferr(:,2),ferr(:,3),'x')
        colorbar
        grid on
            % Create xlabel
            %title('f_3','FontWeight','bold','FontSize',14)
            xlabel({'f_1'},'FontWeight','bold','FontSize',14);
            ylabel({'f_2'},'FontWeight','bold','FontSize',14);
            zlabel({'f_3'},'FontWeight','bold','FontSize',14);
            % Create colorbar
            colorbar('peer',axes1,'FontWeight','bold','FontSize',14);

        %title(['F - ' MOptype])
        hold off
        %figure, scatter3(x(:,1),x(:,2),x(:,3),50,x(:,3),'filled');title(['X - ' MOptype])
        
        xlim([mnf(1),mxf(1)])
        ylim([mnf(2),mxf(2)])
        zlim([mnf(3),mxf(3)])
elseif size(f,2)==4
    figure1=figure;
            axes1 = axes('Parent',figure1,'FontWeight','bold','FontSize',14);
            grid('on');
            hold('all');
        %chin=fill3(ft(:,1),ft(:,2),ft(:,3),'b');
        %alpha(chin,.3);
        %plot3(bt(:,1),bt(:,2),bt(:,3),'.')
        scatter3(fok(:,1),fok(:,2),fok(:,3),50,fok(:,4),'filled')
        plot3(ferr(:,1),ferr(:,2),ferr(:,3),'x')
        %scatter3(ferr(:,1),ferr(:,2),ferr(:,3),40,ferr(:,4))
        colorbar
        title(['color scale - f_4'],'FontWeight','bold','FontSize',14)
        % Create xlabel
            xlabel({'f_1'},'FontWeight','bold','FontSize',14);
            ylabel({'f_2'},'FontWeight','bold','FontSize',14);
            zlabel({'f_3'},'FontWeight','bold','FontSize',14);
            % Create colorbar
            colorbar('peer',axes1,'FontWeight','bold','FontSize',14);
        hold off
        %figure, scatter3(x(:,1),x(:,2),x(:,3),50,x(:,4),'filled');title(['X - ' MOptype])
        
        xlim([mnf(1),mxf(1)])
        ylim([mnf(2),mxf(2)])
        zlim([mnf(3),mxf(3)])
elseif size(f,2)==2
    figure
    %plot(ft(:,1),ft(:,2),'.b')
    hold on
    plot(f(:,1),f(:,2),'.b')
end
end



%Tot_Fcount=mean(fcount);
Tot_Fcount=sum(fcount);

nonConv=find(Cnvrg<1);
[pdom]=filtr_2_pdom(f,f);
nonPP=union(nonConv,pdom);
nnpi=length(nonPP);

ndp=length(Cnvrg);
%  MO_Geom_ProbGraf
 
fok=f;
fok(nonPP,:)=[];
[ev,dl,du,Aregion,tr]=UnifmComput(fok,ptsp);

iok=1:ndp;
%m0=find(du<10^-7);
iok(nonPP)=[];iok(tr)=[];
%iok(du<10^-6)=[];
ndpP=length(iok);
%post_proc(f(iok,:),x,Cnvrg(iok),fcount,Tot_Time,MOptype,ptsp);

TTs=Tot_Time;
TFCs=Tot_Fcount;
eveness=ev;
Par_area1=sum([dl;du])/2/ev;
Par_area=Aregion;
nnonP=nnpi;
fprintf('Results: Time, F Count, Num. NonPareto, Evness, Pareto Area \n')
fprintf('%d %d, %d, %d, %d \n',Tot_Time, Tot_Fcount, nnpi, ev, Par_area)

varargout={TTs,TFCs,eveness,Par_area1,Par_area,nnonP,ndpP};