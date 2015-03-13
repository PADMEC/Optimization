

i=0;
name=['MO2_StatisticsM_SD_2D_' num2str(ceil(cputime))];
ax_id=get(gcf,'CurrentAxes');
while ~isempty(ax_id)
    i=i+1;
    %view(gca,[50,50,50])
    saveas(gcf,[name int2str(i) '.png'])
    saveas(gcf,[name int2str(i) '.fig'])
    %hgsave([name int2str(i)])
    close(gcf)
    ax_id=get(gcf,'CurrentAxes');
end
close(gcf)
display('end')
%Open fig -> save -> close
%i=0;
%for i=1:100
    %open(['filename' int2str(i) '.fig'])
    %saveas(gcf,['filename_' int2str(i) '.png'])
    %close(gcf)
%end