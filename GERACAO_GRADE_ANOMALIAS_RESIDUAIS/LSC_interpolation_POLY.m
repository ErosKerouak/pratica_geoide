clc; clear all; close all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rotina para gerar grade usando Colocação por mínimos quadrados usando função logarítmica 3D %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados=load('remove_pontos.txt');
desvio_g=1;     % mGal

% Importando os pontos da grade interpolada de Anomalia de gravidade residual (mGal)
Dados_grade = load('Grade_interpolar.txt');
resolucao_graus=0.05;    % resolução em graus decimais

% Dados de entrada - coeficientes polinomiais

Xa=readmatrix('Xa_coeficientes polinomiais.txt');

% Dados de saída - grade de dados residuais

nome_arquivo_saida='Grade_AG_res_Poly9.tif';

% Tamanho dos blocos de interpolação em graus
d_lat=0.5;
d_long=0.5;

% Offsets em graus para o efeito de borda nas interpoalções
offset=0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dados das estações
Latitudes=Dados(:,1); Longitudes=Dados(:,2); AG_residuais=Dados(:,4); % mGal
Desvios_g=desvio_g*ones(length(Dados),1); % mGal
Latitudes_interp=Dados_grade(:,1); Longitudes_interp=Dados_grade(:,2);

clear Dados_grade Dados

% Resolvendo por blocos de interpolação
blocos_lat=ceil((max(Latitudes_interp)-min(Latitudes_interp))/d_lat);      % quantidade de blocos em lat
blocos_long=ceil((max(Longitudes_interp)-min(Longitudes_interp))/d_long);  % quantidade de blocos em long
k=0;

for i=1:blocos_lat
for j=1:blocos_long

Lat=Latitudes; Long=Longitudes; AG_res=AG_residuais; N=Desvios_g;
Lati=Latitudes_interp; Longi=Longitudes_interp;

min_lat=min(Lati)+(i-1)*d_lat;
max_lat=min(Lati)+(i)*d_lat;
min_long=min(Longi)+(j-1)*d_long;
max_long=min(Longi)+(j)*d_long;

Lat(Lat<min_lat-offset)=NaN;       % latitude mínima
if i==blocos_lat
Lat(Lat>max_lat+offset)=NaN;       % latitude máxima chegando no limite da área
else
Lat(Lat>=max_lat+offset)=NaN;      % latitude máxima
end
index_lat=~isnan(Lat); index_lat=double(index_lat); index_lat(index_lat==0)=NaN; % criando a máscara para exclusão
Long=Long.*index_lat;  % excluindo das longitudes os ponto fora dos limites de latitude
Long(Long<min_long-offset)=NaN;       % longitude mínima
if j==blocos_long
Long(Long>max_long+offset)=NaN;       % longitude máxima chegando no limite da área
else
Long(Long>=max_long+offset)=NaN;      % longitude máxima
end
index_long=~isnan(Long); index_long=double(index_long); index_long(index_long==0)=NaN; % criando a máscara para exclusão
Lat=Lat.*index_long;   % excluindo das latitudes os pontos fora dos limites de longitude

AG_res=AG_res.*index_long;  % Filtrando os dados para a região de interpolação
N=N.*index_long;            % Filtrando os dados para a região de interpolação

AG_res(isnan(AG_res))=[];   % eliminando os NaN
N(isnan(N))=[];             % eliminando os NaN
Lat(isnan(Lat))=[];         % eliminando os NaN
Long(isnan(Long))=[];       % eliminando os NaN

% Dados da grade de interpolação
Lati(Lati<min_lat)=NaN;       % latitude mínima
if i==blocos_lat
Lati(Lati>max_lat)=NaN;       % latitude máxima chegando no limite da área
else
Lati(Lati>=max_lat)=NaN;      % latitude máxima
end
index_lati=~isnan(Lati); index_lati=double(index_lati); index_lati(index_lati==0)=NaN;
Longi=Longi.*index_lati;
Longi(Longi<min_long)=NaN;       % longitude mínima
if j==blocos_long
Longi(Longi>max_long)=NaN;        % longitude máxima chegando no limite da área
else
Longi(Longi>=max_long)=NaN;      % longitude máxima
end
index_longi=~isnan(Longi); index_longi=double(index_longi); index_longi(index_longi==0)=NaN;

Lati=Lati.*index_longi;     % Filtrando os dados para a região de interpolação

Lati(isnan(Lati))=[];       % eliminando os NaN
Longi(isnan(Longi))=[];     % eliminando os NaN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Colocação 3D e gridagem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gouti=LSCPoly(Longi,Lati,Long,Lat,Xa,N,AG_res);

aloca=(k+1):(max(k))+length(Gouti); aloca=aloca';
k=max(aloca);
Resultado(aloca,:)=[Longi Lati Gouti aloca];

end
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exportar em formato tif %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Resultado=sortrows(Resultado);
nlin=length(unique(Latitudes_interp));
ncol=length(unique(Longitudes_interp));

Imagem=reshape(Resultado(:,3),nlin,ncol);

lat_S_grade=min(Latitudes_interp)-(resolucao_graus/2);
lat_N_grade=max(Latitudes_interp)+(resolucao_graus/2);
lon_W_grade=min(Longitudes_interp)-(resolucao_graus/2);
lon_E_grade=max(Longitudes_interp)+(resolucao_graus/2);
Arquivo=georasterref('RasterSize',size(Imagem),'LatitudeLimits',[lat_S_grade,lat_N_grade],'LongitudeLimits',[lon_W_grade,lon_E_grade]);
geotiffwrite(nome_arquivo_saida,(Imagem),Arquivo)

load handel
sound(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions which fit the covariance function parameters to the data
% by Tiago Lima Rodrigues with base in the Jack (2024)'s griddataLSC function (https://www.mathworks.com/matlabcentral/fileexchange/57342-griddatalsc),
% Developed July 16, 2024.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to perform the least squares collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Least squares collocation using polynomial covariance function
function [SolG]=LSCPoly(Xi,Yi,X,Y,Xa,N,G)
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7); a7=Xa(8);
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5);
s2=(X*ones(size(X))'-ones(size(X))*X').^2+(Y*ones(size(Y))'-ones(size(Y))*Y').^2;
s=sqrt(s2);
% Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4+a5*s.^5+a6*s.^6+a7*s.^7+a8*s.^8+a9*s.^9+a10*s.^10+a11*s.^11)*1D-5;
% Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4+a5*s.^5+a6*s.^6+a7*s.^7)*1D-5;
% Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4+a5*s.^5+a6*s.^6)*1D-5;
Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4);
s2i=(Xi*ones(size(X))'-ones(size(Xi))*X').^2+(Yi*ones(size(Y))'-ones(size(Yi))*Y').^2;
si=sqrt(s2i);
% Csz=(a0+a1*si+a2*si.^2+a3*si.^3+a4*si.^4+a5*si.^5+a6*si.^6+a7*si.^7)*1D-5;
% Csz=(a0+a1*si+a2*si.^2+a3*si.^3+a4*si.^4+a5*si.^5+a6*si.^6)*1D-5;
Csz=(a0+a1*si+a2*si.^2+a3*si.^3+a4*si.^4);
LF=(Czz+diag(N))\eye(size(Czz));
SolG=Csz*(LF*G);
end