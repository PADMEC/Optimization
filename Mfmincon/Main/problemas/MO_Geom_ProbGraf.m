%plot Geometric Problem Grafcs


global funinput coninput

optPoints=funinput{1};
Cc=coninput{1};
Rc=coninput{2};

    iok=1:ndp;
    ierotim=nonPP;
    iok(ierotim)=[];
    fok=f(iok,:);
    ferr=f(nonPP,:);

    %figure,

    figure1=figure;
    %    hold on
        %chin=fill3(ft(:,1),ft(:,2),ft(:,3),'b');
        %alpha(chin,.3);
        %plot3(bt(:,1),bt(:,2),bt(:,3),'.')
        

            axes1 = axes('Parent',figure1,'FontWeight','bold','FontSize',14);
            grid('on');
            hold('all');
    hold on,% title(MOptype)
    tc=0:0.1:2*pi;
    xc=Rc*cos(tc)+Cc(:,1);
    yc=Rc*sin(tc)+Cc(:,2);
    plot(xc,yc,'k','LineWidth',2)
    hold on,
    plot(x(iok,1),x(iok,2),'.b','MarkerSize',15)
    plot(x(ierotim,1),x(ierotim,2),'xk','MarkerSize',15)
    plot(optPoints(:,1),optPoints(:,2),'.r','MarkerSize',30)
            xlabel({'x_1'},'FontWeight','bold','FontSize',14);
            ylabel({'x_2'},'FontWeight','bold','FontSize',14);
            
        xlim([-1,1])
        ylim([-.2,1.6])
        %zlim([mnf(3),mxf(3)])