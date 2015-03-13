
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT RESULTS pushbutton
%
% --- Executes on button press in pushbutton3.
function PlotResults(checkbox,loadstep,varchoice)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Variables shared with POST-PROCESS RESULTS pushbutton function
global lnods coord incrtable availflag Sigs ...
       displ varnames reac ...
       nodepresc varvalues NO_IG hyplasprojname
   
   % MODEFIED NO_IG
    NO_IG=1;
    
    ampl_factor = 20;
    if nargin<2; loadstep=0; end
    if nargin<3; varchoice=0; end 
        
    % Last increment by default
    loadstep(loadstep==0)=incrtable(end);
    % Plot Equivalent Stress by default
    varchoice(varchoice==0)=8;

    % Work out increment whose results are to be plotted 
    % (from INCREMENT NUMBER popupmenu)
    incrid = loadstep; % id number of chosen increment
    if NO_IG
        % Work out variable to be plotted  
        varid = varchoice;
    else
    % (from CONTOUR PLOT VARIABLE popupmenu)
        varid  = varchoice - 1;  % id number of chosen variable
    end
%
%
% Plot meshes
% -----------
%
% Clear axes (current display area) before plotting
handles = figure;

%%% Getting file name
%%fileslas=find(hyplasprojname=='/')+1;
%%titlename = hyplasprojname(fileslas(end):end);

% Creating title
title(hyplasprojname,'Interpreter','none')

%set(gcf,'CurrentAxes',handles.axes3);colorbar('Hide');
%cla; xlim('auto');ylim('auto');
% Clear MESSAGES area
%set(handles.text10,'String','');
hold on;
% Initialise some message flags
incrreacflag=0;incrcontoursflag=0;incrdisplflag=0;
%
% Issue warning if there is no plotting selection
%

if (checkbox<1)
    fprintf(['WARNING: Nothing selected for plotting. '...
                                'The relevant boxes of the '...
                                'PLOTTING... panel must be checked ' ...
                                'in order to plot.']);
else
%
% Plot what is being asked
%
    % Initial mesh
    if checkbox(4)
        patch('Faces',lnods,'Vertices',coord.init,'FaceColor','W','EdgeColor',[0.5 0.5 0.5]);
        % Node numbers
        if checkbox(11)
            npoin=length(coord.init(:,1));
            for ipoin=1:npoin
                text(coord.init(ipoin,1),coord.init(ipoin,2),num2str(ipoin));
            end
        end
    end
    %xdata.deform = xdata.init;ydata.deform = ydata.init;
    coord.deform = coord.init;
    %
    
    %
    % Plot results
    %
    if incrid ~=0
        iincr = incrtable(incrid);   % actual analysis increment number
        %------------------------------------------------------------------
        % Deformed mesh
        %
        if checkbox(3)
            coord.deform = coord.init + ampl_factor * [ displ.x(:,iincr) displ.y(:,iincr) ];
%             patch('Faces',lnods,'Vertices',coord.deform,'FaceColor','W');%'EdgeColor','k',
            % Node numbers
            if checkbox(11)
                npoin=length(coord.deform(:,1));
                for ipoin=1:npoin
                    text(coord.deform(ipoin,1),coord.deform(ipoin,2),num2str(ipoin));
                end
            end
        end
        %------------------------------------------------------------------
        % Plot contours
        %
        if varid ~= 0
            
            % contours
            if checkbox(2)
                %varvalues(varvalues(:,varid,iincr)>100,varid,iincr)=101;
                tfid=get(gca,'title');
                strtitle=get(tfid,'string');
                title('')
                title([strtitle ' - ' varnames{varid}],'Interpreter','none')
%                 patch('Faces',lnods,'Vertices',coord.deform,'FaceVertexCData',...
%                   varvalues(:,varid,iincr),'FaceColor','interp');%,'EdgeColor','k'
                patch('Faces',lnods,'Vertices',coord.deform,'FaceVertexCData',...
                    varvalues(:,varid,iincr),'FaceColor','interp','EdgeColor','none');%
                % PLOTING GP ELEMENT STRESS (for CST)
                [nel,nne]=size(lnods);
                XY=zeros(nel,2);
                for ine=1:nne
                XY=XY+ coord.deform(lnods(:,ine),:)/nne;
                end
                scatter(XY(:,1),XY(:,2),10,Sigs(:,end),'filled');
                colorbar('EastOutside');
            end
            %... and incremental contours
            if checkbox(6)
                if iincr==1
                    incrvarvalues=varvalues(:,varid,iincr);
                else
                   found=0;
                   for iincrprev=iincr-1:-1:1
                       if availflag(3,iincrprev); found=1; break; end
                   end
                   if found
                       incrvarvalues=varvalues(:,varid,iincr)-varvalues(:,varid,iincrprev);
                   else
                       incrvarvalues=varvalues(:,varid,iincr);
                   end
                   if (iincrprev~=iincr-1 || ~found); incrcontoursflag=1; end
                end
                
                patch('Faces',lnods,'Vertices',coord.deform,'FaceVertexCData',incrvarvalues,'FaceColor','interp');%,'EdgeColor','k'
                colorbar('EastOutside');
            end
        elseif varid == 0 ...
               && ( checkbox(2) || checkbox(6))
            fprintf(    ['WARNING: No contour plot variable selected. '...
                        'No contour plot is being plotted. ' ...
                        'To plot contours or incremental contours, '...
                        'a contour plot variable must first be selected.'])
        end
        %------------------------------------------------------------------
        % Plot displacements
        %
        if checkbox(5)
            quiver(coord.init(:,1),coord.init(:,2),displ.x(:,iincr),displ.y(:,iincr),'Color',[0 0 1]);
            hold on;
        end
        %... and incremental displacements
        if checkbox(7)
            if iincr==1
                incrdispl.x=displ.x(:,iincr);incrdispl.y=displ.y(:,iincr);
            else
                found=0;
                for iincrprev=iincr-1:-1:1
                    if availflag(1,iincrprev); found=1; break; end
                end
                if found
                    incrdispl.x=displ.x(:,iincr)-displ.x(:,iincrprev);
                    incrdispl.y=displ.y(:,iincr)-displ.y(:,iincrprev);
                else
                    incrdispl.x=displ.x(:,iincr);
                    incrdispl.y=displ.y(:,iincr);
                end
                if (iincrprev~=iincr-1 || ~found); incrdisplflag = 1; end
            end
            quiver(coord.init(:,1),coord.init(:,2),incrdispl.x,incrdispl.y,'Color',[0 0.5 1]);
            hold on;
        end
        %------------------------------------------------------------------
        % Plot reactions
        %
        if checkbox(1)
            quiver(coord.deform(nodepresc(:,iincr),1),coord.deform(nodepresc(:,iincr),2),reac.x(:,iincr),reac.y(:,iincr),'Color',[1 0 0]);
            hold on;
        end
        %... and incremental reactions
        if checkbox(8)
            if iincr==1
                incrreac.x=reac.x(:,iincr);incrreac.y=reac.y(:,iincr);
            else
                found=0;
                for iincrprev=iincr-1:-1:1
                    if availflag(2,iincrprev); found=1; break; end
                end
                if found
                    incrreac.x=reac.x(:,iincr)-reac.x(:,iincrprev);
                    incrreac.y=reac.y(:,iincr)-reac.y(:,iincrprev);
                else
                    incrreac.x=reac.x(:,iincr);incrreac.y=reac.y(:,iincr);
                end
                if (iincrprev~=iincr-1 || ~found); incrreacflag=1; end
            end
            quiver(coord.deform(nodepresc(:,iincr),1),coord.deform(nodepresc(:,iincr),2),incrreac.x,incrreac.y,'Color',[1 0.5 0]);
            hold on;
        end
        %------------------------------------------------------------------
    else
        fprintf('WARNING: No increments selected for plotting.')
    end
    % Set plot sizes
    xmax=max(xlim);xmin=min(xlim);ymax=max(ylim);ymin=min(ylim);
    lx = xmax-xmin; ly = ymax-ymin;
    lplot = max([lx;ly]);
    xlim([xmin-lplot*0.18 xmin+lplot*1.18]);ylim([ymin-lplot*0.05 ymin+lplot*1.05])
    %
    % Issue mesages about incremental variables
    if incrreacflag && ~incrdisplflag && ~incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental reactions shown are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH REACTIONS ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'reactions for every step and re-run HYPLAS.'])
    elseif ~incrreacflag && incrdisplflag && ~incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental displacements shown are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH DISPLACEMENTS ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'displacements for every step and re-run HYPLAS.'])
    elseif ~incrreacflag && ~incrdisplflag && incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental contours shown are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH CONTOUR DATA ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'variables at nodes for every step and re-run HYPLAS.'])
    elseif incrreacflag && incrdisplflag && ~incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental reactions/displacements shown '...
            'are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH REACTIONS/DISPLACEMENTS DATA ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'reactions/displacements for every step and re-run HYPLAS.'])
    elseif incrreacflag && ~incrdisplflag && incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental reactions/contours shown '...
            'are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH REACTIONS/CONTOUR DATA ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'reactions/variables at nodes for every step '...
            'and re-run HYPLAS.'])
    elseif ~incrreacflag && incrdisplflag && incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental displacements/contours shown '...
            'are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH DISPLACEMENTS/CONTOUR DATA ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'displacements/variables at nodes for every step '...
            'and re-run HYPLAS.'])
     elseif incrreacflag && incrdisplflag && incrcontoursflag
        
        fprintf(...
           ['WARNING: The incremental reactions/displacements/contours '...
            'shown are NOT the '...
            'increment between the current and previous '...
            'load step of the ANALYSIS. They '...
            'are the increment between the last load '...
            'step FOR WHICH REACTION/DISPLACEMENT/CONTOUR DATA '...
            'ARE KNOWN '...
            'and the current step. To plot '...
            'the actual incremental '...
            'values at any step, set your data file to output '...
            'reactions/displacements/variables at nodes for every step '...
            'and re-run HYPLAS.'])       
    end
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%
%%%%%%%%%%%%              END OF functions in               %%%%%%%%%%%%%%%
%%%%%%%%%%%%          RESULTS DISPLAY CONTROL panel         %%%%%%%%%%%%%%%
%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

