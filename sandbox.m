close all

tal = linspace(0, 1, 1000);
N = 7;

R = [(1-tal).^6;
     6*tal.*(1-tal).^5;
     15*tal.^2.*(1-tal).^4;
     20*tal.^3.*(1-tal).^3;
     15*tal.^4.*(1-tal).^2;
     6*tal.^5.*(1-tal);
     tal.^6];

lineStyles = linspecer(N,'qualitative'); 
%lineStyles = linspecer(N,'sequential');
axes('NextPlot','replacechildren', 'ColorOrder',lineStyles);


plot(tal, R(1, :), ...
     tal, R(2, :), ...
     tal, R(3, :), ...
     tal, R(4, :), ...
     tal, R(5, :), ...
     tal, R(6, :), ...
     tal, R(7, :))

legend('R_{0,6}', 'R_{1,6}', 'R_{2,6}', 'R_{3,6}', 'R_{4,6}', 'R_{5,6}', 'R_{6,6}')

xlim([0, 1])
ylim([0, 1])
grid on
set(gcf,'color','w');

C = [1, 2;
     2, 4;
     3, 3;
     4, 1;
     5, 0;
     6, 1;
     7, 2];
 
bezier = C'*R;

figure
plot(bezier(1,:), bezier(2,:), C(:,1), C(:,2), 'r*', C(:,1), C(:,2), 'r*-')
legend('Curva de Bézier', 'Pontos de Controle')
grid on
xlabel('x');
ylabel('y');
set(gcf,'color','w');

C = [1, 2, 0;
     2, 4, 0;
     3, 3, -1;
     4, 1, -1;
     5, 0, -1;
     6, 1, -1;
     7, 2, 0];
 
bezier = C'*R;

figure
plot3(bezier(1,:), bezier(2,:), bezier(3,:))
hold on 
plot3(C(:,1), C(:,2), C(:,3), 'r*')
plot3(C(:,1), C(:,2), C(:,3), 'r*-')
legend('Curva de Bézier', 'Pontos de Controle')
grid on
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gcf,'color','w');


    
