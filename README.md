# inla-r
O INLA (Integrated Nested Laplace Approximations) é um método para inferência bayesiana aproximada que tem sido amplamente utilizado em estatística espacial, especialmente em modelos gaussianos latentes (LGM).

Estatística Espacial:
Envolve a análise de dados que estão associados a localizações geográficas, como mapas de doenças, distribuição de espécies, ou padrões de uso da terra. 
A análise espacial pode incluir a modelagem da dependência espacial entre as observações, ou seja, como a localização de um ponto influencia os valores em outros pontos próximos. 

INLA:
É uma alternativa computacionalmente eficiente para métodos de Monte Carlo de cadeia de Markov (MCMC) em inferência bayesiana, especialmente para modelos com estruturas complexas, como modelos espaciais.
O INLA utiliza uma abordagem chamada de "aproximações de Laplace aninhadas" para calcular distribuições marginais posteriores de parâmetros de interesse.
É implementado no pacote R-INLA do R, uma linguagem e ambiente para computação estatística e gráficos. 

Relação entre Estatística Espacial e INLA:
O INLA tem sido uma ferramenta poderosa para análise estatística espacial, permitindo modelar a dependência espacial de forma eficiente e precisa. 
Ele permite a modelagem de processos espaciais complexos, como aqueles que envolvem equações diferenciais parciais estocásticas (SPDE). 
O INLA tem sido aplicado em diversas áreas da estatística espacial, como mapeamento de doenças, análise de dados de contagem espacial, e estudos de distribuição de espécies. 

Exemplos de aplicação do INLA em Estatística Espacial:

Mapeamento de doenças:
Modelagem de padrões de doenças em áreas geográficas para identificar áreas de maior risco e entender os fatores que influenciam a ocorrência da doença. 

Distribuição de espécies:
Modelagem da distribuição de espécies em um ambiente, considerando fatores ambientais e a interação entre espécies. 

Análise de dados de contagem:
Modelagem de dados de contagem espacial, como o número de ocorrências de um evento em diferentes áreas, levando em conta a influência da localização. 

Bayesian Inference for Multivariate Spatial Models with R-INLA
December 2022
DOI: 10.48550/arXiv.2212.10976
https://www.researchgate.net/publication/366497486_Bayesian_Inference_for_Multivariate_Spatial_Models_with_R-INLA

Pacote "inlabru": https://inlabru-org.github.io/inlabru/

<img width="1200" height="925" alt="Rplot2_" src="https://github.com/user-attachments/assets/f77ed1a5-5575-494d-9163-9a6505093149" />
