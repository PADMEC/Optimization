function ke = keltr3D(area,comp,els,cosn)
    
    ke = zeros(6,6);
    kno = cosn*cosn';
    
    ke = [kno -kno
          -kno kno];
    
    kl = area*els/comp;
    ke=ke*kl;
    
end