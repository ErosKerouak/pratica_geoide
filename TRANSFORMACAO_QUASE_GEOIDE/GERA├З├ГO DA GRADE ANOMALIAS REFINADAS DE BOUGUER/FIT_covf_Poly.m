clc; clear all; close all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Rotina para gerar função covariância e seus parâmetros %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados=load('Pontos Anomalias de gravidade refinada Bouguer.txt');

% grau do polinômio

n=8;

% Dados de saída - coeficientes polinomiais

nome_arquivo_saida='Xa_coeficientes polinomiais';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Função covariância %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lat=Dados(:,1); Long=Dados(:,2); AG_Bouguer=Dados(:,3); clear Dados

[covariance,covdist]=spatialCov(Long,Lat,AG_Bouguer);
covariance=covariance(1:end-1);
covdist=covdist(1:end-1);
s=covdist'; % m

% Ajustamento pelo MMQ - covariance = a0 + a1*s + a2*s^2 + a3*s^3 + a4*s^4 + a5*s^5 + a6*s^6 - melhor estimativa com grau 11
Lb=covariance;  % mGal^2
P=eye(length(covariance),length(covariance));
if n==2
A=[ones(length(Lb),1) s];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2);
covariance_ajust = a0 + a1*s;
elseif n==3
A=[ones(length(Lb),1) s s.^2];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3);
covariance_ajust = a0 + a1*s + a2*s.^2;
elseif n==4
A=[ones(length(Lb),1) s s.^2  s.^3];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3;
elseif n==5
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4;
elseif n==6
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5;
elseif n==7
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6;
elseif n==8
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7); a7=Xa(8);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7;
elseif n==9
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8;
elseif n==10
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9;
elseif n==11
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10;
elseif n==12
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10   s.^11];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10 + a11*s.^11;
elseif n==13
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10   s.^11   s.^12];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);  a12=Xa(13);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10 + a11*s.^11 + a12*s.^12;
elseif n==14
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10   s.^11   s.^12   s.^13];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);  a12=Xa(13);  a13=Xa(14);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10 + a11*s.^11 + a12*s.^12 + a13*s.^13;
elseif n==15
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10   s.^11   s.^12   s.^13   s.^14];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);  a12=Xa(13);  a13=Xa(14);  a14=Xa(15);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10 + a11*s.^11 + a12*s.^12 + a13*s.^13 + a14*s.^14;
elseif n==16
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10   s.^11   s.^12   s.^13   s.^14   s.^15];
Xa=inv(A'*P*A)*A'*P*Lb;
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);  a12=Xa(13);  a13=Xa(14);  a14=Xa(15);  a15=Xa(16);
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10 + a11*s.^11 + a12*s.^12 + a13*s.^13 + a14*s.^14 + a15*s.^15;
end

R2=1 - sum((covariance - covariance_ajust).^2)/sum((covariance - mean(covariance)).^2)

% s=s(1:end-1);
% covariance_ajust=covariance_ajust(1:end-1);
% covariance=covariance(1:end-1);

figure(1);
plot(s,covariance_ajust,'DisplayName','Fitted covariance function - polynomial','LineWidth',1.5);   hold on
plot(s,covariance,'DisplayName','Empirical covariance function','LineWidth',1.5); hold off
str=['R^2 = ' num2str(R2)];    dim = [.6 .45 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'BackgroundColor','white');
grid on; grid minor;  lgd = legend; lgd.NumColumns = 1;
xlabel('Distance (°)','FontName','Times New Roman','FontSize',14)
ylabel('Covariance (mGal^2)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exportação dos resultados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gerando txt para exportar
Xa_table = array2table(Xa);
writetable(Xa_table,nome_arquivo_saida,'Delimiter',' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions which fit the covariance function parameters to the data
% by Jack (2024). griddataLSC (https://www.mathworks.com/matlabcentral/fileexchange/57342-griddatalsc),
% MATLAB Central File Exchange. Retrieved January 13, 2024.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to find the 2-D Empirical Covariance of the input data
function [covariance,covdist]=spatialCov(X,Y,G)
% smax=sqrt((max(X)-min(X))^2+(max(Y)-min(Y))^2); % maximum range of the data
smax=4; % maximum range of the data
ds=sqrt((2*pi*(smax/2)^2)/(length(X)));         % average spatial density
covariance=zeros(round(smax/ds)+2,1);
ncov=zeros(round(smax/ds)+2,1);
for i = 1:length(G)
        for j = 1:length(G)
            if j~=i
          xx=(X(i)-X(j));
          yy=(Y(i)-Y(j)); 
          r = sqrt(xx.^2+yy.^2);
          ir = round((r/ds))+1;
          if r<smax
            covariance(ir) = covariance(ir) + G(i)*G(j) ;
            ncov(ir) = ncov(ir)+1;
          end
            end
        end
end
for i = 1:length(covariance)
  if ncov(i)~=0
      covariance(i) = covariance(i)/ncov(i);
  end
end  
covdist=(0:length(covariance)-1)*ds;
end