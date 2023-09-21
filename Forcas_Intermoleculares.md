
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

Vamos considerar primeiramente um exemplo extremamente simples de dois átomos interagentes de Hélio e calcular $E_{\rm int}$ em algumas distâncias interatômicas distintas $R$. Você usará o pacote gratuito Psi4 para calcular as energias total que você precisa para realizar a subtração. Quando fizermos estes calculos para uma série de valores de $R$, você será capaz de construir as famosas *curva de energia potencial* - o gráfico de $E_{\rm int}(R)$ em função de $R$.

BELEZA, mas como você escolherá o método de estrutura eletrônica para calcular $E_{\rm A-B}$, $E_{\rm A}$, and $E_{\rm B}$? Vamos começar pela escolha mais simples e tentar o método Hartree-Fock (HF). Caso o HF não seja acurado o suficiente, tentaremos um método mais sofisticado como o *coupled-cluster* com excitações simples, duplas e triplas perturbativas - CCSD(T). Se você nunca ouviu falar sobre CCSD(T) antes, vamos apenas dizer que **(1)** geralmente é muito preciso (muitos se referem a este método como *gold standard* da teoria de estrutura eletrônica) e **(2)** é computacionalmente muito pesado para moléculas maiores - de certa forma até inviável. Com relação ao conjunto de funções de base, vamos escolher a base de Dunning aug-cc-pVTZ que deve estar OK para HF e CCSD(T).

## Dímero de Hélio

Aqui está um script python que emprega o Psi4 para computar a curva de energia potencial para dois átomos de Hélio empregando uma abordagem supramolecular. Primeiramente vamos importar algumas bibliotecas importantes
``` 
import time
import numpy as np
import scipy
from scipy.optimize import *
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)
import psi4
import matplotlib.pyplot as plt

# Set Psi4 & NumPy Memory Options
psi4.set_memory('2 GB')
psi4.core.set_output_file('output.dat', False)

numpy_memory = 2

psi4.set_options({'basis': 'aug-cc-pVTZ',
              'e_convergence': 1e-10,
              'd_convergence': 1e-10,
              'INTS_TOLERANCE': 1e-15})

```
