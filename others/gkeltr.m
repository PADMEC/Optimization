	function ke = gkeltr(area,ls,els,teta,kelem,jelem,ielem)
	ke = zeros(4,4);

s = sin(teta);
c = cos(teta);
kl = els/ls;
            
if kelem == ielem
kt=[kl*c*c kl*c*s
    kl*c*s kl*s*s];

ke=[kt -kt
    -kt kt];

end;                  
