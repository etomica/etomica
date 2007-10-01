package etomica.potential;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.box.Box;
import etomica.math.SpecialFunctions;
import etomica.nbr.CriterionNone;
import etomica.nbr.NeighborCriterion;
import etomica.space.IVector;
import etomica.space3d.IVector3D;

public class EwaldSummation {
	
	
	/*
	 * Ewald Sum 
	 * 
	 * This class is developed based on Ewald Summation to efficiently sum the 
	 * interaction between an ion and all its periodic images
	 * 
	 * Refereces:
	 * 1. Allen and Tildesley, Computer Simulation of Liquids, 1st Ed., Pg 156-160
	 * 2. U. Essmann, L. Perera, and M. L. Berkowitz, J. Chem. Phys. 103(19), 1995 
	 * 		[specifically from Equations 2.1 - 2.5]
	 * 
	 * @param
	 * 
	 * nVectorMax : is number of the vector of the periodic images
	 * alpha : is the coefficient that determine the contributions of each terms in the sum
	 * 		   usually is set to 5/L, where L is the length of the simulation box
	 * 
	 * @author: Tai Tan
	 */
	
	
	public EwaldSummation(Box box, AtomAgentManager atomAgentManager, int nVectorMax, double alpha){
		
		this.box = box;
		this.atomAgentManager = atomAgentManager;
		this.moleculeBasis = new AtomSetSinglet();
		this.nVectorMax = nVectorMax;
		this.alpha = alpha;
		atomPair = new AtomPair();
		
		setCriterion(new CriterionNone());
		int rCut = nVectorMax*nVectorMax;
		
		int numNVector = 0;
		
		for (int i=-nVectorMax; i<nVectorMax; i++){
			for (int j=-nVectorMax; j<nVectorMax; j++){
				for (int k=-nVectorMax; k<nVectorMax; k++){
					
					int sq = i*i + j*j + k*k;
					if (sq !=0 && sq <= rCut){
						
						numNVector += 1;
					}
				}
			}
		}
		
		nVector = new IVector3D[numNVector];

		numNVector = 0;
		
		for (int i=-nVectorMax; i<nVectorMax; i++){
			for (int j=-nVectorMax; j<nVectorMax; j++){
				for (int k=-nVectorMax; k<nVectorMax; k++){
					
					int sq = i*i + j*j + k*k;
					if (sq !=0 && sq <= rCut){
						
						nVector[numNVector] = (IVector3D)box.getSpace().makeVector(new double[] {i,j,k});
						numNVector += 1;
					}
				}
			}
		}
		
	}
	
	/*
	 *  Ewald Summation for the Real Part
	 *  
	 */
	public double EwaldSumReal(){
		
		int numNVector = nVector.length;
		double uReal = 0.0;
		
		AtomSet atomList = box.getLeafList();
		int numAtom = atomList.getAtomCount();
		IVector rij = box.getSpace().makeVector();
		IVector rijNv = box.getSpace().makeVector();
		
		for (int vecCounter=-1; vecCounter<numNVector; vecCounter++){
			
			for (int i=0; i<numAtom; i++){
				IAtom atomi = atomList.getAtom(i);
				IVector posAtomi = ((IAtomPositioned)atomi).getPosition();
				double chargei = ((MyCharge)atomAgentManager.getAgent(atomi)).charge;
				atomPair.atom0 = atomi;
				
				
				for (int j=0; j<numAtom; j++){
					IAtom atomj = atomList.getAtom(j);
					atomPair.atom1 = atomj;
					
					if((atomi == atomj || criterion.accept(atomPair))&& i==-1){
						continue;
					}
					
					IVector posAtomj = ((IAtomPositioned)atomj).getPosition();
					double chargej = ((MyCharge)atomAgentManager.getAgent(atomj)).charge;
					
					rij.Ev1Mv2(posAtomj, posAtomi);
					
					if (vecCounter == -1){
						rijNv.E(rij);
					} else {
						rijNv.Ev1Pv2(rij, nVector[vecCounter]);
					}
					double rijNvMagnitude = Math.sqrt(rijNv.squared());
					
					uReal += 0.5*chargei*chargej*SpecialFunctions.erfc(alpha*rijNvMagnitude)/ rijNvMagnitude;
				}
			}
		}
		
		return uReal;
	}
	
	/*
	 *  Ewald Summation for the Fourier Part
	 *  
	 */
	public double EwaldSumFourier(){

		int numNVector = nVector.length;
		
		/*
		 * Computing uFourier
		 */
		
		double preFactor = 1/(2*Math.PI*box.volume());
		double L = Math.pow(box.volume(), 1/3)  ;
		double b = Math.PI*Math.PI/((alpha*alpha)*(L*L));
		double sumVectorTerm = 0.0;
		
		for (int vecCounter=0; vecCounter<numNVector; vecCounter++){
			
			// L^2*[exp(-Pi^2*n^2/ (L^2*alpha^2))/ n^2]
			double nSquared = nVector[vecCounter].squared();
			double Q = L*L*Math.exp(-b*nSquared)/nSquared;
			
			
			/*
			 * Solve expression for S(n)*S(-n)
			 */
			
			AtomSet atomList = box.getLeafList();
			int numAtom = atomList.getAtomCount();
			double pl = 2*Math.PI/L;
			
			// S(n) and S(-n)
			
			double realMagnitude = 0.0;
			double imagMagnitude = 0.0;
			
			for (int i=0; i<numAtom; i++){
				IAtom atomi = atomList.getAtom(i);
				IVector posAtomi = ((IAtomPositioned)atomi).getPosition();
				double chargei = ((MyCharge)atomAgentManager.getAgent(atomi)).charge;
				
				double Sn = pl*nVector[vecCounter].dot(posAtomi);
				
				realMagnitude += chargei* Math.cos(Sn);
				imagMagnitude += chargei* Math.sin(Sn);
			}	
				
				/*
				 *  S(n)*S(-n)
				 *  They are conjugated to each other
				 *  
				 *  for (a+ib)*(a-ib) = a^2 + b^2
				 */
			sumVectorTerm += Q*(realMagnitude*realMagnitude + imagMagnitude*imagMagnitude);
		}
		
		return preFactor*sumVectorTerm;
		
	}
	
	
	/*
	 *  Allen and Tildesley, Computer Simulation of Liquids, 1st Ed.(1987)
	 *   Pg 159, Eq 5.20
	 */
	
	public double EwaldSumSelf(){
		// Get the molecule list of the first species
		// What you really want is all the lists of the molecules
		AtomSet moleculeList = ((IAtomGroup)box.getSpeciesMaster().getAgentList().getAtom(0)).getChildList();
		//
		
		/*
		 * molecules can be monoatomic or multiatomic
		 */
		int numMolecules = moleculeList.getAtomCount();
		
		double prefcancel = alpha/Math.sqrt(Math.PI);
		double cancelTerm = 0.0;
		double pairTerm = 0.0;
		
		
		for (int i=0; i<numMolecules; i++){
			IAtom molecule = moleculeList.getAtom(i);
			
			if (!(molecule instanceof IAtomGroup)){
				/*
				 * Monatomic
				 */
				double chargei = ((MyCharge)atomAgentManager.getAgent(molecule)).charge;
				cancelTerm += prefcancel*chargei*chargei;
			}
		 
		
			else {
				
				/*
				 * Multi atomic
				 */
				
				IVector dabVector = box.getSpace().makeVector();
				int numSites = ((IAtomGroup)molecule).getChildList().getAtomCount();
			
				// cancel-Term
				for (int a=0; a<numSites; a++){
					IAtom sitea = ((IAtomGroup)molecule).getChildList().getAtom(a);
					
					double chargea = ((MyCharge)atomAgentManager.getAgent(sitea)).charge;
					cancelTerm += prefcancel*chargea*chargea;
				}
				
				// self-Term
				moleculeBasis.atom = molecule;
				iterator.setBasis(moleculeBasis);
				iterator.reset();
				
				for (AtomSet pair = iterator.next(); pair!= null; pair = iterator.next()){
					IAtom sitea = pair.getAtom(0);
					IAtom siteb = pair.getAtom(1);
					
					IVector posAtoma = ((IAtomPositioned)sitea).getPosition();
					IVector posAtomb = ((IAtomPositioned)siteb).getPosition();
					
					double chargea = ((MyCharge)atomAgentManager.getAgent(sitea)).charge;
					double chargeb = ((MyCharge)atomAgentManager.getAgent(siteb)).charge;
				
					dabVector.Ev1Mv2(posAtomb, posAtoma);
					double dab = Math.sqrt(dabVector.squared());
					
					pairTerm += chargea*chargeb*(1-SpecialFunctions.erfc(alpha*dab))/dab;
					}
				
				}
		
			}
		
		return cancelTerm + 0.5*pairTerm;
	}
	
	public double EwaldSum(){
		
		return EwaldSumReal() + EwaldSumFourier() + EwaldSumSelf();
	}
	
	public void setCriterion(NeighborCriterion criterion){
		this.criterion = criterion;
	}
	
	public void setBondedIterator(AtomsetIteratorBasisDependent iterator){
		this.iterator = iterator;
		
	}
	
	protected final AtomAgentManager atomAgentManager;
	protected final IVector3D[] nVector;
	protected final int nVectorMax;
	protected final double alpha;
	protected final Box box;
	protected NeighborCriterion criterion;
	protected final AtomPair atomPair;
	protected AtomsetIteratorBasisDependent iterator;
	protected final AtomSetSinglet moleculeBasis;
	
	
	public static class MyCharge{
		
		public MyCharge(double charge){
			this.charge = charge;
		}
		
		public final double charge;
		
	}
}
