	function ke = keltr(area,ls,els,teta)
	ke = zeros(4,4);

cs = [cos(teta), sin(teta)];
kt = cs'*cs;

kl = area*els/ls;

ke=kl*[kt -kt
    -kt kt];
