function PlasticMove(Fs,Us,Sigs,type,aviname)

 %%%%%%%%%%%% VIDEO PLASTIFCTN %%%%%%%%%%%
    nsteps=size(Us,2);
    im=0;clear('MOV');
    for i=1:nsteps
        im=im+1;
        plot_sig_vet(Fs(:,i),Us(:,i),Sigs(:,i),type);
        MOV(im)=getframe(gcf);
        close(gcf)
    end
    movie(MOV,3,1)
        MOVIE2AVI(MOV,aviname)
    if nargin>4
        MOVIE2AVI(MOV,aviname,'COMPRESSION','Indeo5')
    end