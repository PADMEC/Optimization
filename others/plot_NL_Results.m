function id_plot=plot_NL_Results(Fs,Us,Sigs,str_title,id_plot)


if nargin<4
    str_title='';
end

    ms=max(abs(Sigs));
    mu=max(abs(Us));

    if nargin>4
        figure(id_plot{1});
        hold on, h=plot(mu,ms,'.-');
    else
        id_plot{1}=figure;
        h=plot(mu,ms,'r.-');
    end
    xlabel('desloc')
    ylabel('sig')
    set(h,'DisplayName',str_title);
    legend('show')
    %legend(h,str_title);
    
    iF=find(Fs(:,1),1);
    F1=abs(Fs(iF,:));
    

    if nargin>4
        figure(id_plot{2});
        hold on, h=plot(F1,mu,'.-');
    else
        id_plot{2}=figure;
        h=plot(F1,mu,'r.-');
    end
    %hold on, plot(Fs(iF,:),ms/max(ms),'r.-')
    xlabel('F')
    ylabel('desloc')
    legend(h,str_title);
    
    if nargin>4
        figure(id_plot{3});
        hold on, h=plot(F1,ms,'.-');
    else
        id_plot{3}=figure;
        h=plot(F1,ms,'r.-');
    end
    xlabel('F')
    ylabel('sig')
    legend(h,str_title);

    nsteps=size(Us,2);
    for i=1:nsteps
        nU(i)=norm(Us(:,i));
    end
    if nargin>4
        figure(id_plot{4});
        hold on, h=plot(F1,nU,'.-');
    else
        id_plot{4}=figure;
        h=plot(F1,nU,'r.-');
    end
    xlabel('F')
    ylabel('norm(U)')
    legend(h,str_title);
        