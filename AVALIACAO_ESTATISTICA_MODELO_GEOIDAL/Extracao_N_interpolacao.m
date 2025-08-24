clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades
Grade=double(imread('D:\AULAS MÉTODOS FÍSICOS\TRABALHO PRÁTICO 2025-1\RESTORE\Modelo_geoidal.tif'));

resolucao_graus=0.05;      % Resolução espacial normal em °

latN=-23.475-(resolucao_graus/2);     % latN
lonW=-53.525+(resolucao_graus/2);    % longW

Pontos_RN_GPS=readmatrix('Pontos_RN_GPS.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extração por interpolação %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lin,col]=size(Grade);

latS=(latN-(lin-1)*resolucao_graus);
lonE=(lonW+(col-1)*resolucao_graus);

grade_lat=(latN:-resolucao_graus:latS)'.*ones(1,col);
grade_long=(lonW:resolucao_graus:lonE).*ones(lin,1);

latP=Pontos_RN_GPS(:,1); longP=Pontos_RN_GPS(:,2);

NP=interp2(grade_long,grade_lat,Grade,longP,latP,'linear');