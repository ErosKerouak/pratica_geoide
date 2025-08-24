clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades
Grade=double(imread('MDS_MERIT_SRTM15PLUS_900m_fill.tif'));

resolucao_graus=0.0083333333;      % Resolução espacial do MDS em °

latN=-21.2504166667-(resolucao_graus/2);     % latN
lonW=-55.7495833333+(resolucao_graus/2);    % longW

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

HP=interp2(grade_long,grade_lat,Grade,longP,latP,'linear');