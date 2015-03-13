%Generete domained non-domined GRAPHCS

%Connect all Pareto points in FT, before
 FT=fok;
% FT=[FT;fok];

[pdom,pndom]=filtr_2_pdom(fok,FT);

[length(pdom) length(pndom)]

figure,scatter3(fok(pndom,1),fok(pndom,2),fok(pndom,3),50,fok(pndom,3),'filled')
hold on,plot3(fok(pdom,1),fok(pdom,2),fok(pdom,3),'xr')
xlim([150, 250]),ylim([2000, 6000]),zlim([0, 0.25])
view([90,90,90])
xlabel(cc(1,:)),ylabel(cc(2,:)),zlabel(cc(3,:)),title(metext)
