%get_gauss

    availflag(3,iincr)=1; % set contour plot variables availability flag
    fgetl(fid); fgetl(fid); fgetl(fid);
    ipoincounter=0;tline='*';breakflag=0;
    while ~( ipoincounter == npoin && breakflag==1 )
        tline = fgetl(fid);
        if tline==-1;
            varvalues(ipoin,:,iincr) = nodvarval;
            break;
        end
        if length(tline)>11 && strcmp(tline(1:12),' Node number')
            firsttimeflag=firsttimeflag+1;
            ipoin=str2num(tline(13:21));ipoincounter=ipoincounter+1;
            nodvarval=[];
        elseif  strcmp(tline,'')
            %varvalues(ipoin,:,iincr) = nodvarval;
            if ipoincounter==npoin
                breakflag=1;
            end
        elseif ipoincounter==1 && firsttimeflag==1 % get variable names only 
                                                   % in the first time the
                                                   % relevant data is read
            pointers = findstr('=',tline);
            a=sscanf(tline,'%*s = %g ');
            nodvarval= [ nodvarval a'];
            for ii=1:length(pointers)
                varnames = [ varnames tline(pointers(ii)-6:pointers(ii)-1)];
                %nodvarval= [ nodvarval str2num(tline(pointers(ii)+2:pointers(ii)+13))];
            end
        else
            a=sscanf(tline,'%*s = %g ');
            if ~isstr(a)
                nodvarval= [ nodvarval a'];
            end
%                 pointers = findstr('=',tline);
%                 for ii=1:length(pointers)
%                     nodvarval= [ nodvarval str2num(tline(pointers(ii)+2:pointers(ii)+13))];
%                 end
        end
    end
