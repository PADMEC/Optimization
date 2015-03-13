function write_Zfile(POD_Z,filename)


fid = fopen(filename, 'w');
%Writing binnarilly!
fwrite(fid, size(POD_Z),'uint');fwrite(fid, POD_Z,'double');
fclose(fid);

%[nl,nc]=size(POD_Z);
% fprintf(fid,'%d  %d \n',[nl,nc]);
% format='\n';
% for i=1:nc
%     format=['%11.9g ' format];
% end
% fprintf(fid,format,POD_Z');