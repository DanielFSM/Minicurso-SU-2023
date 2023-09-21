
## Introdução 
Esta atividade de laboratório foi pensada para auxiliar os estudantes matriculados no minicurso "Brincando com Forças Intermoleculares" da Universidade de Brasília. Este minicurso consiste em uma série de laboratórios computacionais sobre forças intermoleculares fracas e o calculo/interpretação da energia de interação entre duas moléculas.
A energia de interação pode ser decomposta em contribuições fisicamente significativas (eletrostática, indução, dispersão e troca) usando a teoria da perturbação adaptada à simetria (SAPT). Neste exercício, calcularemos as energias de interação completas e sua decomposição SAPT usando os procedimentos do pacote de software Psi4, processando e analisando os dados com NumPy e Matplotlib.

Objetivos de Aprendizagem:
1. Reconhecer e apreciar a ubiquidade e diversidade das interações intermoleculares.
2. Comparar e contrastar os métodos supramoleculares e perturbativos para calcular as energias de interação.
3. Analisar e interpretar as contribuições SAPT na forma eletrostática, indução dispersão e troca em diferentes separações intermoleculares.

# Forças Intermoleculares Fracas 

Nesta atividade, você irpa examinar as algumas proprieades das interações fracas entre moléculas. Como as subunidades moleculares não estão conectados por nenhuma ligação covalente ou iônica, nos referimos usualmente as interações intermoleculares como in *interações não-covalentes*. Suponha que queremos calcular a energia de interação entre a molécula A e a molécula B para uma certa geometria do complexo A-B (óbviamente, esta energia de interação depende de quão distante as moléculas estão e como elas estão orientadas). A forma mais simples de realizar esta tarefa seria através da subtração (esta abordadem é chamada de *abordagem supramolecular*, é a mais empregada por químicos teóricos):


$E_{\rm int}=E_{\rm A-B}-E_{\rm A}-E_{\rm B}$

em que $E_{\rm X}$ é a energia total do sistema X, computado usando nossa teoria favorita de estrutura eletrônica e um certo conjunto de funções de base. Um valor negativo de  $E_{\rm int}$ indica que a molécula A e B tem energia mais baixa quando estão juntas naquela configuração do que quando elas estão infinitamente separadas, então eles formam um complexo ligado que pode ser estável pelo menos em temperaturas muito baixas. Um valor positivo de $E_{\rm int}$ indica que o complexo A-B não é ligado - é energeticamente mais favorável se A e B se separarem (repulsões). 

Let's consider a simple example of two interacting helium atoms and calculate $E_{\rm int}$ at a few different interatomic distances $R$. You will use Psi4 to calculate the total energies that you need to perform subtraction. When you do so for a couple different $R$, you will be able to sketch the *potential energy curve* - the graph of $E_{\rm int}(R)$ as a function of $R$.

OK, but how should you pick the electronic structure method to calculate $E_{\rm A-B}$, $E_{\rm A}$, and $E_{\rm B}$? Let's start with the simplest choice and try out the Hartree-Fock (HF) method. In case HF is not accurate enough, we will also try the coupled-cluster method with single, double, and perturbative triple excitations - CCSD(T). If you haven't heard about CCSD(T) before, let's just state that it is **(1)** usually very accurate (it's even called the *gold standard* of electronic structure theory) and **(2)** very expensive for larger molecules. For the basis set, let's pick the augmented correlation consistent triple-zeta (aug-cc-pVTZ) basis of Dunning which should be quite OK for both HF and CCSD(T).
