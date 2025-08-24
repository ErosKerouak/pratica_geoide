clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades
Grade=double(imread('.\GERAÇÃO DA GRADE ANOMALIAS REFINADAS DE BOUGUER\Grade_AG_Bouguer_Poly12.tif'));

resolucao_graus=0.05;      % Resolução espacial da grade de AG_Bouguer em °

latN=-22.475-(resolucao_graus/2);     % latN
lonW=-54.525+(resolucao_graus/2);    % longW

Pontos_grade=readmatrix('Pontos_calculo_grade.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extração por interpolação %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lin,col]=size(Grade);

latS=(latN-(lin-1)*resolucao_graus);
lonE=(lonW+(col-1)*resolucao_graus);

grade_lat=(latN:-resolucao_graus:latS)'.*ones(1,col);
grade_long=(lonW:resolucao_graus:lonE).*ones(lin,1);

latP=Pontos_grade(:,1); longP=Pontos_grade(:,2);

AG_BouguerP=interp2(grade_long,grade_lat,Grade,longP,latP,'linear');