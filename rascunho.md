

\chapter*{Resumo}
\addcontentsline{toc}{chapter}{Resumo}



---

# Prompt — Escrita Científica Clara e Eficaz

> Você é um assistente de escrita científica. Sua missão é reescrever, revisar ou sugerir melhorias em textos acadêmicos seguindo estritamente as ordens abaixo:
>
> * **Comunique, não apenas apresente.** Garanta que o leitor entenda a ideia, não apenas veja dados.
> * **Mantenha sujeito e verbo juntos.** Nunca esconda o verbo atrás de interrupções longas.
> * **Coloque a informação importante no fim da frase.** A posição final (“stress position”) deve carregar o que é novo ou crucial.
> * **Comece cada frase com o já conhecido.** Use a posição inicial (“topic position”) para dar contexto e continuidade.
> * **Dê a cada sentença uma função única.** Uma frase = um ponto; um parágrafo = uma ideia central.
> * **Expresse a ação no verbo.** Nunca enterre a ação em substantivos abstratos.
> * **Forneça contexto antes da novidade.** Prepare o leitor antes de apresentar conceitos ou resultados novos.
> * **Alinhe a estrutura com a importância.** A forma deve refletir o peso relativo das ideias.
> * **Explicite todas as conexões lógicas.** Nunca obrigue o leitor a adivinhar relações entre frases ou parágrafos.
> * **Obedeça às expectativas antes de quebrá-las.** Só viole princípios se souber exatamente o efeito pretendido.
>
> Sua saída deve ser clara, direta e estruturada. Ao revisar, mostre **primeiro o texto corrigido** e depois, se necessário, uma breve justificativa com base nas ordens acima.

---

Com base nos capítulos abaixo escreva o resumo.
O Resumo deve determinar CLARA e SUCINTAMENTE qual o assunto e quais os principais resultados do trabalho.

-O Resumo deve conter (i) os principais objetivos do trabalho, (ii) a metodologia empregada, (iii) os resultados mais importantes e (iv) as conclusões mais relevantes;
-Comece definindo o problema investigado e a metodologia empregada em poucas frases, então selecione os resultados e as conclusões ais úteis do trabalho.
- O Resumo NÃO deve ter citação de Figuras, Tabelas ou Referências Bibliográficas.

```LaTex

\chapter{Introdução}

O campo de gravidade da Terra fornece a base física para transformar alturas geométricas em alturas físicas e, com isso, integrar posicionamento GNSS ao nivelamento clássico. Essa integração depende de modelos confiáveis do geoide e do quase-geoide, cuja acurácia regional condiciona aplicações em engenharia, cartografia e monitoramento geodinâmico. Em escala continental, modelos satelitais do geopotencial descrevem bem os comprimentos de onda longos; em escala local, os dados terrestres resolvem feições curtas e heterogêneas. O desafio surge na junção dessas escalas: combinar medições irregulares, ruidosas e afetadas pela topografia com a informação espectral dos modelos globais, de modo estatística e fisicamente consistente.

Para enfrentar esse desafio, adota-se uma estratégia em duas frentes e um ponto de encontro. Na frente espectral, remove-se do sinal observado a contribuição de um modelo global do geopotencial, reduzindo a variância de \emph{background} e isolando as componentes residuais relevantes à área de estudo. Na frente espacial, modela-se a dependência entre observações por uma função de covariância adequada e aplica-se Colocação por Mínimos Quadrados (LSC) para interpolar o sinal residual em grade regular, controlando ruído e bordas. O ponto de encontro ocorre na integração geodésica: a função de Stokes, na forma modificada compatível com a remoção espectral, converte anomalias residuais de gravidade em ondulação residual do geoide, à qual se agrega a contribuição indireta da topografia (efeito indireto) para obter superfícies fisicamente interpretáveis.

Neste trabalho, segue-se essa arquitetura de solução com escolhas explícitas e verificáveis: (i) definem-se anomalias residuais de Helmert a partir de dados terrestres, após subtração da contribuição de um modelo geopotencial de referência (GOCO06s, até grau $300$) e do efeito indireto de Helmert; (ii) estima-se a covariância empírica como função da separação angular e ajusta-se um polinômio isotrópico de grau $9$ até $4^\circ$, produzindo um modelo suave e operacional para uso em LSC; (iii) realiza-se a interpolação por LSC com modelagem explícita do ruído de observação e blocagem espacial para mitigar condicionamento e efeitos de borda; (iv) integra-se a grade residual por PVCG com função de Stokes modificada (truncagem compatível com o grau removido), com tratamento da singularidade central e aproximação analítica da zona interna; (v) calcula-se o efeito indireto em pontos e em grade a partir de MDE de alta resolução e, por fim, (vi) extraem-se valores de geoide e quase-geoide em marcos de rede de nivelamento (RNs) para validação por resíduos, viés, desvio-padrão e RMSE.

A partir desse desenho, formulam-se as seguintes \textbf{hipóteses de trabalho}: (H1) um modelo de covariância isotrópico polinomial até $4^\circ$ captura a estrutura espacial necessária à LSC para o campo residual considerado; (H2) a truncagem espectral da função de Stokes, compatível com o grau removido do modelo global, produz uma ondulação residual coerente com os comprimentos de onda resolvidos pelos dados terrestres; (H3) a modelagem de ruído ao nível de $\sim 1{,}0$\,mGal estabiliza a solução sem introduzir suavização excessiva; (H4) a inclusão do efeito indireto da topografia reduz enviesamentos sistemáticos na comparação com RNs.

Com base nessas hipóteses, enuncia-se o \textbf{problema} que orienta o estudo: \emph{como obter, a partir de anomalias de Helmert em pontos irregularmente distribuídos e de um modelo geopotencial de referência, uma superfície de ondulação geoidal/quase-geoidal residual que seja estatística e fisicamente consistente com marcos de rede de nivelamento na área de estudo, explicitando o papel do modelo de covariância e da modificação de Stokes?}

O \textbf{objetivo geral} é construir e validar um modelo residual de geoide/quase-geoide para a área de estudo combinando dados terrestres e informação satelital de referência. Como \textbf{objetivos específicos}, pretende-se: (i) padronizar e calcular anomalias residuais de gravidade; (ii) ajustar uma função de covariância polinomial e avaliar sua aderência; (iii) interpolar o campo residual por LSC com controle de ruído e bordas; (iv) integrar por Stokes modificado para obter a ondulação residual; (v) estimar o efeito indireto e expressá-lo nas grandezas de interesse; e (vi) validar o quase-geoide perante RNs por análise de resíduos e métricas-síntese.

\chapter{Discussão}


O conjunto \texttt{remove\_pontos.csv} foi aceito para uso subsequente porque passou por controle de consistência com erro ponto a ponto e métricas (erro absoluto máximo e RMSE) calculadas a partir da comparação com a planilha de referência, conforme relatado em \emph{Resultados}; as figuras \Figauto{fig:hist-dg-res} e \Figauto{fig:mapa-dg-res} sustentam essa decisão ao sintetizarem, respectivamente, a amplitude global e a cobertura espacial de $\Delta g_{\text{res}}$, condições necessárias para estimar de forma estável a covariância empírica e suportar o ajuste polinomial posterior.


O ajuste de grau $n=9$ reproduziu $C_{\text{emp}}(s)$ no intervalo $s \le 4^\circ$ com $R^{2}=0{,}9549$ (\Figauto{fig:covariance-fit}), indicador objetivo de aderência entre modelo e dados; esse desempenho justifica, com base no número reportado, o uso dessa função na construção de $\mathbf{C}_{zz}$ e $\mathbf{C}_{sz}$ para a colocação, mantendo a hipótese operacional de isotropia e a escala de correlação implícita no polinômio ajustado.


A estratégia de LSC adotou $\Sigma_n=\mathrm{diag}(1{,}0~\mathrm{mGal})$ e blocagem com margem, escolhas explicitadas em \emph{Resultados} para controlar condicionamento e bordas; o produto gerado (\texttt{Grade\_AG\_res\_Poly9.tif}, $0{,}05^\circ$) cumpre o objetivo de prover um campo contínuo de entrada para a integração, e sua construção se apoia diretamente no ajuste polinomial validado por $R^2$, de modo que a estabilidade do interpolador se encontra fundamentada no desempenho observado da covariância ajustada.


A grade $N_{\text{res}}$ resultou da integração de $\Delta g_{\text{res}}$ com domínio recuado de $1^\circ$, remoção explícita da singularidade e aproximação analítica da zona interna, decisões metodológicas registradas em \emph{Resultados}; o mapa temático, o histograma e o perfil (\Figauto{fig:nres-mapa}–\Figauto{fig:nres-perfil}) documentam, sem ambiguidade, a continuidade espacial do campo, a distribuição de valores e a variação ao longo de uma seção latitudinal, fornecendo a base empírica para qualquer avaliação subsequente de amplitude e gradientes do modelo residual.


Os cálculos pontuais de EI e a grade $EI_N$ (\texttt{EI\_N.tif}) foram produzidos com raio de integração de $100$~km, correção de curvatura e interpolação linear com \emph{nearest} como \emph{fallback}, parâmetros explicitados em \emph{Resultados}; as figuras \Figauto{fig:ei-hist}, \Figauto{fig:ei-pontos} e \Figauto{fig:ei-grid} mostram, respectivamente, a amplitude, a variabilidade espacial nos pontos de cálculo e o campo na malha do modelo, o que permite, com base nesses dados, avaliar a contribuição topográfica indireta tanto local quanto regional antes de qualquer combinação com $N_{\text{res}}$.


A validação do quase-geoide em $n=32$ RNs apresentou viés de $-0{,}434$~m, desvio-padrão de $0{,}144$~m e RMSE de $0{,}456$~m, números reportados em \emph{Resultados} e que caracterizam, de forma objetiva, o desempenho do modelo frente às observações; por definição estatística, a remoção do viés reduziria o RMSE ao nível do desvio-padrão (aprox.\ $0{,}144$~m), o que quantifica, apenas com base nos dados, o ganho potencial de uma correção média simples. As figuras \Figauto{fig:zeta-hist}, \Figauto{fig:zeta-mapa}, \Figauto{fig:zeta-scatter} e \Figauto{fig:zeta-trends} sustentam essas afirmações ao exibirem a distribuição dos resíduos, sua organização espacial, a relação 1:1 entre $\zeta_{\text{model}}$ e $\zeta_{\text{RN}}$ e possíveis dependências com latitude/longitude; adicionalmente, a \Figauto{fig:nmodel-rn} documenta o campo $N_{\text{model}}$ nos RNs usados, permitindo cotejar, nos dados, cobertura e desempenho local.


Considerados em conjunto, os resultados mostram (i) dados pontuais de $\Delta g_{\text{res}}$ com cobertura e amplitude documentadas por histograma e mapa, (ii) função de covariância ajustada com alta aderência ($R^{2}=0{,}9549$) que fundamenta a LSC, (iii) grade interpolada de $\Delta g_{\text{res}}$ apta à integração, (iv) modelo $N_{\text{res}}$ caracterizado por mapa, histograma e perfil, e (v) EI quantificado e espacializado, além de (vi) validação do quase-geoide com viés e dispersão explicitamente mensurados. Essas constatações derivam diretamente das métricas e figuras apresentadas e, portanto, suportam decisões práticas como a aplicação de correção de viés na validação, testes de sensibilidade ao ruído $\Sigma_n$ e ao grau do polinômio, e o uso combinado de $N_{\text{res}}$ e $EI_N$ em etapas subsequentes; quaisquer interpretações geofísicas finas (por exemplo, causas de padrões regionais específicos) devem ser baseadas nas cartas e estatísticas já incluídas, evitando extrapolações além do que os dados mostram.

\chapter{Conclusão}

Este trabalho demonstra um encadeamento eficaz entre remoção espectral, modelagem estocástica e integração geodésica para construir um modelo residual de geoide a partir de anomalias de Helmert e de um modelo geopotencial de referência. Parte-se do problema de compatibilizar escalas espectrais distintas -- comprimentos de onda longos bem descritos por modelos globais e feições curtas resolvidas por dados terrestres -- e chega-se a superfícies residuais fisicamente interpretáveis e mensuradas contra marcos de referência.

No nível estocástico, confirma-se que um modelo de covariância polinomial isotrópico até $4^\circ$ captura a estrutura espacial relevante do campo residual: o ajuste de grau $9$ alcança $R^{2}=0{,}9549$, indicando aderência suficiente para sustentar a Colocação por Mínimos Quadrados. No nível numérico, a LSC, com modelagem explícita de ruído ao nível de $1{,}0\,\mathrm{mGal}$ e blocagem com margem, produz uma grade contínua e estável para integração, sem sinais de suavização excessiva nos diagnósticos apresentados. No nível geodésico, a integração por PVCG com função de Stokes modificada -- truncada de modo compatível com o grau removido -- gera uma ondulação residual contínua e coerente com os comprimentos de onda preservados, com tratamento explícito da singularidade e aproximação analítica da zona interna.

A avaliação externa por marcos de rede de nivelamento mostra desempenho mensurável do quase-geoide modelado: viés de $-0{,}434$\,m, desvio-padrão de $0{,}144$\,m e $\mathrm{RMSE}=0{,}456$\,m para o conjunto analisado. Por definição, a simples remoção do viés reduziria o $\mathrm{RMSE}$ ao patamar do desvio-padrão ($\sim 0{,}144$\,m), quantificando de forma objetiva o ganho potencial de uma correção média. A estimativa e a espacialização do efeito indireto (EI) ficam estabelecidas como componente adicional do fluxo; embora os resultados indiquem contribuição topográfica sistemática, a combinação explícita $N_{\text{res}}+EI_N$ na validação final permanece uma etapa a consolidar para testar integralmente a hipótese de redução de enviesamento.

À luz das hipóteses formuladas, conclui-se que: (H1) é suportada pelos resultados do ajuste e pela estabilidade da LSC; (H2) é corroborada pelos produtos integrados e pelos diagnósticos de continuidade; (H3) é consistente com o equilíbrio observado entre condicionamento e preservação de variabilidade; (H4) permanece parcialmente aberta, exigindo validação com a soma explícita do EI ao campo residual na comparação com RNs.

As principais limitações decorrem de escolhas operacionais que impõem compromissos entre resolução e estabilidade: isotropia do modelo de covariância, janela angular de $4^\circ$, grau polinomial fixado, truncagem espectral alinhada ao modelo global, e número/modo de distribuição dos RNs disponíveis. Adicionalmente, aproximações geométricas (área de pixel e correção de curvatura) e a incerteza do MDE afetam a estimativa do EI em curtas escalas.
```
