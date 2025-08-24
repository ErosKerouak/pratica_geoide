clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades
Ag_res=double(imread('Grade_AG_res_Poly9.tif'))*1D-5;         % anomalias de gravidade residuais em m2/s2 com resolução espacial normal

resolucao_graus=0.05;      % Resolução espacial normal em °
raio_integracao_graus=1;       % Raio de integração em °
grau_modificacao=300;

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
latN=-22.47500000000000141-(resolucao_graus/2);     % latN
lonW=-54.5249999999999986+(resolucao_graus/2);    % longW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grades de coordenadas geodésicas da área de integração para a grade de anomalias residuais e MDE
[lin,col]=size(Ag_res);

latS=(latN-(lin-1)*resolucao_graus);
lonE=(lonW+(col-1)*resolucao_graus);

grade_lat=(latN:-resolucao_graus:latS)'.*ones(1,col)*pi/180; % rad - grade de latitudes da anomalia residual
grade_long=(lonW:resolucao_graus:lonE).*ones(lin,1)*pi/180;  % rad - grade de longitudes da anomalia residual

% Pontos de integração
latQ=grade_lat(:);       longQ=grade_long(:);      AgQ=Ag_res(:);
sin_latQ=sin(latQ);     cos_latQ=cos(latQ);

% Pontos de cálculo
latP_N=latN-raio_integracao_graus;   % latitude norte da grade de cálculo do modelo, descontando o raio de integração
lonP_W=lonW+raio_integracao_graus;   % longitude oeste da grade de cálculo do modelo, descontando o raio de integração
latP_S=latS+raio_integracao_graus;   % latitude sul da grade de cálculo do modelo, descontando o raio de integração
lonP_E=lonE-raio_integracao_graus;   % longitude leste da grade de cálculo do modelo, descontando o raio de integração

lin_modelo=floor((latP_N-latP_S)/resolucao_graus)+1;  % número de colunas da grade do modelo - arredondado para baixo para ser inteiro, +1 para pegar o último valor do intervalo
col_modelo=floor((lonP_E-lonP_W)/resolucao_graus)+1;  % número de linhas da grade do modelo - arredondado para baixo para ser inteiro, +1 para pegar o último valor do intervalo

grade_latP=(latP_N:-resolucao_graus:latP_S)'.*ones(1,col_modelo)*pi/180; % rad - grade de latitudes da anomalia residual
grade_longP=(lonP_W:resolucao_graus:lonP_E).*ones(lin_modelo,1)*pi/180;  % rad - grade de longitudes da anomalia residual

latP=grade_latP(:);       longP=grade_longP(:);

% Anomalia de gravidade e altitude normal dos pontos de cálculo P
AgP=interp2(grade_long,grade_lat,Ag_res,longP,latP,'cubic');

% Gravidade normal no Elipsoide dos pontos de cálculo P
sin_latP=sin(latP);
gama0=(9.7803267715*(1+0.0052790414*(sin_latP.^2)+0.0000232718*(sin_latP.^4)+0.0000001262*(sin_latP.^6)+0.0000000007*(sin_latP.^8)));  % m/s2

% Constantes
res_rad=resolucao_graus*pi/180;              % Resolução em rad
raio_integ_rad=raio_integracao_graus*pi/180; % raio de integração em rad
R=6371000;                                   % Raio da esfera de mesmo volume do Elpsoide em metros
res_rad_div_2=(res_rad/2);
cte=4*pi*R*gama0;
n_mod=2:grau_modificacao;
cte_mod=(2*n_mod+1)./(n_mod-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Solução do PVCG Helmert-Stokes %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_res=zeros(length(latP),1); % alocando espaço para ganhar velocidade

% Áreas
Ak=abs((R^2)*res_rad*(sin(latQ+res_rad_div_2)-sin(latQ-res_rad_div_2)));

% Contribuição da zona mais interna
Ni=(R./gama0).*sqrt((cos(latP)*res_rad^2)/pi).*AgP; % m

for p=1:length(latP)

% Filtrando os pontos de integração dentro da área de integração
cos_psi=sin_latQ*sin_latP(p)+cos_latQ*cos(latP(p)).*cos(longQ-longP(p));
psi=real(acos(cos_psi));

psi(psi>raio_integ_rad)=NaN;  % valores de psi maiores que o raio de integração em rad são transformados em NaN
index=~isnan(psi); % cria um vetor de índice do mesmo tamanho do vetor psi onde na posição de psi igual a NaN index é igual a zero, o contrário igual a 1
index=double(index); % passando de logical para double (número)
index(index==0)=NaN; % onde tiver 0 no index transforma em NaN

% filtrando todos os vetores:
psi(isnan(psi))=[]; % retirando do vetor psi os NaN
cos_psi_integracao=cos_psi.*index; cos_psi_integracao(isnan(cos_psi_integracao))=[];
latQ_integracao=latQ.*index; latQ_integracao(isnan(latQ_integracao))=[];
longQ_integracao=longQ.*index; longQ_integracao(isnan(longQ_integracao))=[];
AgQ_integracao=AgQ.*index; AgQ_integracao(isnan(AgQ_integracao))=[];
Ak_integracao=Ak.*index; Ak_integracao(isnan(Ak_integracao))=[];

% Detectando a localização do ponto de cálculo nos vetores de coordenadas
L1=find(abs((latP(p)-latQ_integracao))==min(abs((latP(p)-latQ_integracao))));
L2=find(abs((longP(p)-longQ_integracao))==min(abs((longP(p)-longQ_integracao))));
pos_p=intersect(L1,L2);

% Função de Stokes com modificação de Wang e Gore (1969)

s2=sin((latP(p)-latQ_integracao)./2).*sin((latP(p)-latQ_integracao)./2)+sin((longP(p)-longQ_integracao)./2).*sin((longP(p)-longQ_integracao)./2).*cos(latP(p)).*cos(latQ_integracao);
s=s2.^0.5;
S_psi=1./s-4-6*s+10*s2-(3-6*s2).*log(s2+s);
S_psi(pos_p)=NaN;   % desconsiderando S-psi do ponto de cálculo (infinito)

Pn = legendreN(grau_modificacao,cos_psi_integracao); % Polinômios de Legendre de grau máximo n e x = cos(psi)
S_psi_SH=cte_mod*Pn(:,3:end)';                       % começando de 3 para dispensar os termos de grau 0 e 1

% Altura geoidal residual em metros
N_res(p)=sum(((S_psi-S_psi_SH').*Ak_integracao.*AgQ_integracao)/cte(p),'omitnan')+Ni(p);

end

toc

Resultado=[latP*180/pi longP*180/pi N_res ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Função rápida dos Polinômios de Legendre de grau n e ordem zero %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pn = legendreN(nmax,t)
Pn = zeros(length(t), nmax+1);
Pn(:, 1) = 1;    
Pn(:, 2) = t; 
for n = 2:nmax 
    Pn(:, n+1) = (1 - n) / n * Pn(:, n - 1) + (2 * n - 1) / n * t .* Pn(:, n);
end
end