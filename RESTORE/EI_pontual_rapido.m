clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Importando o MDE
H = double(imread('MDS_MERIT_SRTM15PLUS_900m_fill.tif'));
[lin,col]=size(H);

% Coordenadas geodésicas das estações de cálculo
Pontos=load('Pontos_calculo_grade.txt');

% Resolução
res=0.008333333333333333868;   % Resolução do MDE em °
raio_integracao=100000;        % raio de integração em metros

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
latN=-21.2504166666666663-(res/2);    % latN
lonW=-55.7495833333333337+(res/2);   % longW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grade de coordenadas geodésicas do MDE
gradeMDA_lat=(latN:-res:(latN-(lin-1)*res))'.*ones(1,col);
gradeMDA_long=(lonW:res:(lonW+(col-1)*res)).*ones(lin,1);

% Coordenadas geodésicas, altitude ortométrica e gravidade normal no elipsoide das estações de cálculo
latP=Pontos(:,1); longP=Pontos(:,2);

HP=interp2(gradeMDA_long,gradeMDA_lat,H,longP,latP,'cubic');    % interpolação para obtenção do valor de HN do MDE nas posições de cálculo

gama0=(9.7803267715*(1+0.0052790414*(sin(latP*pi/180).^2)+0.0000232718*(sin(latP*pi/180).^4)+0.0000001262*(sin(latP*pi/180).^6)+0.0000000007*(sin(latP*pi/180).^8)));

% Constantes
G = 6.67259*1E-11;   % m3 kg-1 s-2
rho=2670;          % kg m-3
res_metros=res*108000;
EI=zeros(length(latP),1);
cte0=-pi*G*rho*HP.^2;    cte1=-(G*rho/6)*res_metros^2;    cte2=(3*G*rho/40)*res_metros^2;
metros=3600*30; % graus para metros
denominador=3*6371000^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Efeito indireto da Topografia no potencial %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(latP)    % percorrendo os n pontos de cálculo

lo=(((gradeMDA_lat-latP(ii,1))*metros).^2+((gradeMDA_long-longP(ii,1))*metros).^2).^0.5; % distância plana em metros

lo=lo(:); lo(lo>raio_integracao)=NaN; index=~isnan(lo); index=double(index); index(index==0)=NaN; % filtrando a distância pelo raio e criando índice de máscara
lo(isnan(lo))=[];             % eliminando os NaN
lo=lo+(lo.^3)/denominador;    % correção de curvatura na distância

H_prisma=H(:).*index; H_prisma(isnan(H_prisma))=[]; % aplicando a máscara de raio de integração e eliminando os NaN

EI(ii,:)=cte1*sum((H_prisma.^3-HP(ii,1)^3)./lo.^3) + cte2*sum((H_prisma.^5-HP(ii,1).^5)./lo.^5); % EI no potencial em m2/s2
end

EI_N=EI./gama0  % Efeito indireto em N em m
EI_Ag=(EI+cte0)*(0.3086)./gama0; % Efeito indireto em AG em mGal

toc