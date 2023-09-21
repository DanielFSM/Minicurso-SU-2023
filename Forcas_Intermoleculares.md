
## Forças Intermoleculares 
This lab activity is designed to teach students about weak intermolecular interactions, and the calculation and interpretation of the interaction energy between two molecules. The interaction energy can be broken down into physically meaningful contributions (electrostatics, induction, dispersion, and exchange) using symmetry-adapted perturbation theory (SAPT). In this exercise, we will calculate complete interaction energies and their SAPT decomposition using the procedures from the Psi4 software package, processing and analyzing the data with NumPy and Matplotlib.

Prerequisite knowledge: the Hartree-Fock method, molecular orbitals, electron correlation and the MP2 theory. The lab also assumes all the standard Python prerequisites of all Psi4Education labs.

Learning Objectives: 
1. Recognize and appreciate the ubiquity and diversity of intermolecular interactions.
2. Compare and contrast the supermolecular and perturbative methods of calculating interaction energy.
3. Analyze and interpret the electrostatic, induction, dispersion, and exchange SAPT contributions at different intermolecular separations.

Author: Konrad Patkowski, Auburn University (patkowsk@auburn.edu; ORCID: 0000-0002-4468-207X)

Copyright: Psi4Education Project, 2020

# Weak intermolecular interactions 

In this activity, you will examine some properties of weak interactions between molecules. As the molecular subunits are not connected by any covalent (or ionic) bonds, we often use the term *noncovalent interactions*. Suppose we want to calculate the interaction energy between molecule A and molecule B for a certain geometry of the A-B complex (obviously, this interaction energy depends on how far apart the molecules are and how they are oriented). The simplest way of doing so is by subtraction (in the so-called *supermolecular approach*):

\begin{equation}
E_{\rm int}=E_{\rm A-B}-E_{\rm A}-E_{\rm B}
\end{equation}

where $E_{\rm X}$ is the total energy of system X, computed using our favorite electronic structure theory and basis set. A negative value of $E_{\rm int}$ means that A and B have a lower energy when they are together than when they are apart, so they do form a weakly bound complex that might be stable at least at very low temperatures. A positive value of $E_{\rm int}$ means that the A-B complex is unbound - it is energetically favorable for A and B to go their separate ways. 

Let's consider a simple example of two interacting helium atoms and calculate $E_{\rm int}$ at a few different interatomic distances $R$. You will use Psi4 to calculate the total energies that you need to perform subtraction. When you do so for a couple different $R$, you will be able to sketch the *potential energy curve* - the graph of $E_{\rm int}(R)$ as a function of $R$.

OK, but how should you pick the electronic structure method to calculate $E_{\rm A-B}$, $E_{\rm A}$, and $E_{\rm B}$? Let's start with the simplest choice and try out the Hartree-Fock (HF) method. In case HF is not accurate enough, we will also try the coupled-cluster method with single, double, and perturbative triple excitations - CCSD(T). If you haven't heard about CCSD(T) before, let's just state that it is **(1)** usually very accurate (it's even called the *gold standard* of electronic structure theory) and **(2)** very expensive for larger molecules. For the basis set, let's pick the augmented correlation consistent triple-zeta (aug-cc-pVTZ) basis of Dunning which should be quite OK for both HF and CCSD(T).