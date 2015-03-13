function Pareto_Plot(fok,ferr,Objectivesi,MOptype)
    
global legnd
%f=[f(:,3) f(:,1) f(:,2) ];

    f=[fok;ferr];

    mxf=(max(f)).*1.1;
    mnf=(min(f)).*.9;
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
                scatter3(fok(:,1),fok(:,2),fok(:,3),50,fok(:,1),'filled')
                plot3(ferr(:,1),ferr(:,2),ferr(:,3),'x')
                colorbar
                grid on
                % Create Labels
                if length(Objectivesi)==3
                    xlabel(Objectivesi{1},'FontWeight','bold','FontSize',14);
                    ylabel(Objectivesi{2},'FontWeight','bold','FontSize',14);
                    zlabel(Objectivesi{3},'FontWeight','bold','FontSize',14);
                else
                    xlabel({'f_1'},'FontWeight','bold','FontSize',14);
                    ylabel({'f_2'},'FontWeight','bold','FontSize',14);
                    zlabel({'f_3'},'FontWeight','bold','FontSize',14);
                end
                % Create colorbar
                colorbar('peer',axes1,'FontWeight','bold','FontSize',14);

            title(MOptype,'FontWeight','bold','FontSize',14)
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
        MOptTypes={'ws','mmx','nbi','nnc','nbim','nncm'};
        method=find(strcmpi(MOptype,MOptTypes));
        sizep=0;newfig=0;
        if isempty(method),method=-1;end
        % Set options
        switch method
            case 1, pt='*b';legnd={'WS'};newfig=1;
            case 2, pt='^g';legnd{end+1}='Min Max';
            case 3, pt='ok';legnd{end+1}='NBI';sizep=1;
            case 4, pt='+r';legnd{end+1}='NNC';sizep=1;
            case 5, return;%pt='.k';leg{i}='NBIm';
            case 6, return;%pt='.r';leg{i}='NNCm';
            otherwise,newfig=1;sizep=2;pt='.';
                legend(legnd{:})
        end
        % Open new figure
        if newfig
            fig1=figure;
            axes('Parent',fig1,'FontWeight','bold','FontSize',14);
            hold on                
        end
        
        % Set Labls
        if length(Objectivesi)==2
            xlabel(Objectivesi{1},'FontWeight','bold','FontSize',14);
            ylabel(Objectivesi{2},'FontWeight','bold','FontSize',14);
        else
            xlabel({'f_1'},'FontWeight','bold','FontSize',14);
            ylabel({'f_2'},'FontWeight','bold','FontSize',14);
        end
        % Set point sizes
        if sizep==0
            sizep=[8,10,12];
        elseif sizep==1
            sizep=[8,10,12];
        elseif sizep==2;
            sizep=[15,8,12];
        end
        plot(f(:,1),f(:,2),pt,'MarkerSize',sizep(1),'LineWidth',2)
%         plot(fok(:,1),fok(:,2),'.','MarkerSize',sizep(2),'LineWidth',2)
%         if ~isempty(fok),legnd{end+1}='Pareto';end
        plot(ferr(:,1),ferr(:,2),['xk'],'MarkerSize',sizep(3),'LineWidth',2)
        if ~isempty(ferr),legnd{end+1}='nonPareto';end
    end