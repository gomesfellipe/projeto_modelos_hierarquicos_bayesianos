# Modelos lineares hierárquicos bayesianos

**Identificação do Projeto**: Projeto de Iniciação à Pesquisa

# Introdução

Modelar um fenômeno aleatório consiste em realizar aﬁrmações sobre o processo
gerador dele. Muitas vezes essas aﬁrmações são realizadas em vários níveis como, por
exemplo, quando as observaçõoes pertencem a grupos diferentes e cada grupo tem suas
próprias propriedades (média, variância, entre outras). Nesses casos recorre-se a modelos
hierárquicos, que também são conhecidos como modelos multiníveis.

Aplicaço˜es desses modelos podem ser encontradas em várias áreas tais como na
Educaçao, nas Ciências Sociais e na Saúde. Suponha que demógrafos desejam examinar
como diferenças no desenvolvimento da economia nacional podem interferir na relação
entre o grau educacional dos adultos e a taxa de fertilidade. Para isso, pode-se
utilizar 2 estágios: nível nacional (indicadores econˆomicos) e nível domiciliar (educaçao e
fertilidade). Ou suponha que o interesse está em medir o rendimento escolar dos alunos e,
para isso, utiliza-se 4 estágios: os alunos, as turmas, as escolas e os órgao administradores
ou a regiao.

Modelos de regressão linear explicam a variável resposta através de variáveis
explicativas e supo˜e que dada as variáveis explicativas, as variáveis respostas são
independentes. Modelos lineares hierárquicos sao generalizaço˜es dos modelos de regressão
linear pois assumem que as observaço˜es das unidades pertencentes ao agregado sao
dependentes.

Ao modelar a variável resposta, há o surgimento de variáveis desconhecidas que sao
chamadas de parâmetros. A inferencia bayesiana permite que qualquer crença que se
tenha sob esse parâmetro seja considerada através da distribuiçao a priori em conjunto
com a função de verossimilhança.

# Objetivo

O objetivo deste trabalho foi o estudo sobre a modelagem hierárquica
bayesiana. Para isso, dados foram simulados a ﬁm de analisar a
sensibilidade da distribuição a priori atribuída e aprender sobre como inferir sobre os
parâmetros desconhecidos.

# Motivação

A inferência bayesiana permite que crenças adquiridas antes da amostragem sejam
consideradas na modelagem do problema. E, quando não há crença, há a possibilidade
de uma distribuiçao nao informativa ser atribúıda tendo resultados similares a inferência
clássica.

Muitas vezes modelar hierarquicamente facilita a atribuição das distribuiço˜es e a
interpretação dos parâmetros

# Metodologia

A inferência sobre os parâmetros desconhecidos será realizada através da abordagem
Bayesiana. Sendo assim, será necessário atribuir uma distribuiçao a priori ao conjunto
de parâmetros desconhecidos para combinar com a função de verossimilhança e, assim,
obter a distribuição a posteriori. A distribuiçao a posteriori pode ser muito ou pouca
afetada pela distribuiçao a priori. Além disso, dependendo das informaço˜es advindas da
amostra, há a necessidade de colocarmos uma distribuição mais informativa para alguns
parâmetros.

Para estimar os parâmetros, recorreremos aos métodos de Monte Carlo via cadeias de
Markov (MCMC). Como a distribuiçao conjunta a posteriori costuma nao possuir forma
analítica conhecida, será utilizado o Amostrador de Gibbs com passos de Metropolis
Hastings para a estimaçao.
Sendo assim, o aluno simulará dados e fará uma análise de sensibilidade para investigar
as questoes já apontadas e aprenderá a utilizar o MCMC.

# Relato

Este projeto teve início por volta de setembro de 2017 e seu principal
objetivo foi o estudo sobre a modelagem hierárquica sob a perspectiva de
inferência bayesiana pois esta metodologia foi adotada posteriormente em
minha monografia, de título "Uma aplicação de modelo hierárquico
bayesiano na modelagem da dor em rescém nascidos submetidos à punção de
calcâneo" e como objetivo secundário, devido à complexidade da
metodologia estudada, foi planejada uma revisão da literatura de modelos
de regressão linear simples sob o paradigma clássico e sob o paradigma
bayesiano para o aprendizado de novos elementos relacionados ao seu
ajuste.

Para isso foram utilizadas referências bibliográficas sugeridas pela
professora Dra. Patrícia para a compreensão da metodologia e toda a
implementação da parte computacional referentes à: simulação,
manipulação e implementação do algorítimo estudado foram realizados com
o uso do software R (versão 1.0.153) que é uma linguagem e um ambiente
para programação estatística e além disso todo o texto produzido foi
feito em \LaTeX, uma linguagem para a produção de artigos científicos.

Toda semana ocorriam reuniões com a orientadora para discutir o tema,
retirar eventuais dúvidas dos cálculos realizados e no decorrer do
projeto quando fez-se nescessário as simulações em R esses encontros
também envolviam a avaliação do andamento da produção dos códigos
computacionais.

Diversos resultados interessantes foram obtidos com o estudo de tanto no
estudo do modelo de regressão linear simples quanto o hieárquico. No
modelo de regressão linear simples além dos dados simulados ainda foi
utilizado um conjunto de dados reais que contam com registros de tempo e
distância até a parada de um carro ao apertar o freio em alta velocidade
e no modelo de regressão linear hierárquico todos os parâmetros
utilizados para simular a amostra foram recuperados.

De maneira resumida o projeto foi muito benéfico não apenas
acrescentando o conhecimento sobre uma nova metodologia para o ajuste de
modelos de regressão linear como também gerou diversos benefícios
secundários como: aperfeiçoar habilidades de simulação dos dados,
treinar a escrita de artigos cientificos com \LaTeX, aprimorar a
habilidade em programação em linguagem R melhorando assim tanto a minha
formação como reconhecendo o verdadeiro valor da educação e da ciência.

# Bibliograﬁa

Banerjee, S., Gelfand, A. E. e Carlin, B. P. (2003) Hierarchical Modeling and Analysis
for Spatial Data. Chapman & Hall/CRC.

Gamerman, D. e Lopes, H. F. (2006) Monte Carlo Markov Chain: Stochastic
Simulation for Bayesian Inference. London: Chapman & Hall, second edn.

Gelman, A. e Hill, J. (2006) Data Analysis Using Regression and Multilevel /
Hierarquical Models. Cambridge University Press.

Raudenbush, S. W., Bryk, A. S. (2001) Hierarchical linear models: applications and
data analysis methods. Sage.

Souza, R. de O. Modelagem do nível de dor e estresse de recém-nascidos internados
em UTI neonatal utilizando um modelo hierárquico Bayesiano com dados longitudinais.
Universidade Federal Fluminense, RJ, Brasil: Trabalho de Conclusao de Curso
Bacharelado em Estatística, 2015.

Velarde, L. G. C., Migon, H. S. , Alcoforado, D. A. Hierarchical Bayesian models
applied to air surveillance radars. European Journal of Operational Research, v. 184, p.
1155-1162, 2008


