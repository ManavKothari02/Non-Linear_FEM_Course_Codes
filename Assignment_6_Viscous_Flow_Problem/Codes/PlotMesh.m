function PLOTMESH(NNM,GLXY)
 S = NNM;
 X = GLXY(1:NNM,1);
 Y = GLXY(1:NNM,2);
 [XX, YY] = meshgrid(X,Y);
 figure;
 for I = 1:S
     LABELS(I)=I;
 end
 LABELS = num2cell(LABELS);
 plot(X,Y,'o');
 text(X,Y,LABELS,'VerticalAlignment','bottom','HorizontalAlignment','right')
 axis equal;
 end