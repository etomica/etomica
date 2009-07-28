package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMolecular;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.MoleculeSetSinglet;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.math.SpecialFunctions;
import etomica.nbr.CriterionNone;
import etomica.nbr.NeighborCriterion;
import etomica.space.Space;

public class EwaldSummation implements IPotentialMolecular {
	
	
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


	public double energy(IMoleculeList atoms) {
		double energy = EwaldSum();
		System.out.println("Energy Ewald Sum: "+ energy);
		return energy;
	}

	public double getRange() {
		// TODO Auto-generated method stub
		return Double.POSITIVE_INFINITY;
	}

	public int nBody() {
		// TODO Auto-generated method stub
		return 0;
	}

	public void setBox(IBox box) {
		// do nothing
		
	}

	public EwaldSummation(IBox box, AtomLeafAgentManager atomAgentManager, int nVectorMax, Space _space){
		
		this.box = box;
		this.space = _space;
		this.atomAgentManager = atomAgentManager;
		this.moleculeBasis = new MoleculeSetSinglet();
		this.nVectorMax = nVectorMax;
		this.alpha = 5/Math.pow(box.getBoundary().volume(), 1.0/3.0);
		atomPair = new AtomPair();
		
		setCriterion(new CriterionNone());
		int rCut = nVectorMax*nVectorMax;
		
		int numNVector = 0;
		
		for (int i=-nVectorMax; i<nVectorMax+1; i++){
			for (int j=-nVectorMax; j<nVectorMax+1; j++){
				for (int k=-nVectorMax; k<nVectorMax+1; k++){
					
					int ii = Math.abs(i);
					int jj = Math.abs(j);
					int kk = Math.abs(k);
					
					if (ii > 0){
						ii--;
					}
					
					if (jj > 0){
						jj--;
					}
					
					if (kk > 0){
						kk--;
					}
					
					int check = ii*ii + jj*jj + kk*jj;
					int sq = i*i + j*j + k*k;
					
					if (sq !=0 && check <= rCut){
						
						numNVector += 1;
					}
				}
			}
		}
		
		nVector = new IVectorMutable[numNVector];

		numNVector = 0;
		
		for (int i=-nVectorMax; i<nVectorMax+1; i++){
			for (int j=-nVectorMax; j<nVectorMax+1; j++){
				for (int k=-nVectorMax; k<nVectorMax+1; k++){
					
					int ii = Math.abs(i);
					int jj = Math.abs(j);
					int kk = Math.abs(k);
					
					if (ii > 0){
						ii--;
					}
					
					if (jj > 0){
						jj--;
					}
					
					if (kk > 0){
						kk--;
					}
					
					int check = ii*ii + jj*jj + kk*jj;
					int sq = i*i + j*j + k*k;
					
					if (sq !=0 && check <= rCut){
						//System.out.println("[i, j, k]: "+i+" " +j+" " + k);
						nVector[numNVector] = space.makeVector(new double[] {i,j,k});
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
		
		IMoleculeList moleculeList = box.getMoleculeList();
		//
		
		/*
		 * molecules can be monoatomic or multiatomic
		 */
		int numMolecules = moleculeList.getMoleculeCount();
		
		int numNVector = nVector.length;
		double uReal = 0.0;
		
		IVectorMutable rij = space.makeVector();
		IVectorMutable rijNv = space.makeVector();
		
		for (int vecCounter=-1; vecCounter<numNVector; vecCounter++){
			
			for (int i=0; i<numMolecules; i++){
				IMolecule moleculei = moleculeList.getMolecule(i);
				
					
					
				/*
				 * Multi atomic for Ewald Real
				 */
				
				int numSites = moleculei.getChildList().getAtomCount();
					
				for (int a=0; a<numSites; a++){
					IAtom sitea = moleculei.getChildList().getAtom(a);
					IVectorMutable posAtoma = ((IAtomPositioned)sitea).getPosition();
					double chargea = ((MyCharge)atomAgentManager.getAgent(sitea)).charge;
					atomPair.atom0 = sitea;
					
					for (int j=0; j<numMolecules; j++){
						IMolecule moleculej = moleculeList.getMolecule(j);
						
						for (int b=0; b<numSites; b++){
							IAtom siteb = moleculej.getChildList().getAtom(b);
							IVectorMutable posAtomb = ((IAtomPositioned)siteb).getPosition();
							double chargeb = ((MyCharge)atomAgentManager.getAgent(siteb)).charge;
							atomPair.atom1 = siteb;
						
							rij.Ev1Mv2(posAtomb, posAtoma);
							
							if (vecCounter == -1){
								if (i==j && (a==b||criterion.accept(atomPair))){
									continue;
								}
								rijNv.E(rij);
							} else {
								rijNv.E(nVector[vecCounter]);
								rijNv.TE(box.getBoundary().getBoxSize());
								rijNv.PE(rij);
							}
							
							double rijNvMagnitude = Math.sqrt(rijNv.squared());
							
							
							uReal += 0.5*chargea*chargeb*SpecialFunctions.erfc(alpha*rijNvMagnitude)/ rijNvMagnitude;
							//System.out.println("uReal difference: " + (uReal-uRealBefore));
							//System.out.println("rijMagnitude: "+ rijNvMagnitude);
						}	
					}
				}
				//System.out.println("Molecule i: "+ moleculei);
				//System.out.println("Ureal: "+ uReal);
				
			} // loop over the molecules
			
		}// nVector Loop
		//System.exit(1);
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
		
		double preFactor = 1/(2*Math.PI*box.getBoundary().volume());
		double L = Math.pow(box.getBoundary().volume(), 1.0/3.0)  ;
		double b = Math.PI*Math.PI/((alpha*alpha)*(L*L));
		double sumVectorTerm = 0.0;
		
		for (int vecCounter=0; vecCounter<numNVector; vecCounter++){
			
			// L^2*[exp(-Pi^2*n^2/ (L^2*alpha^2))/ n^2]
			double nSquared = nVector[vecCounter].squared();
			double Q = L*L*Math.exp(-b*nSquared)/nSquared;
			
			
			/*
			 * Solve expression for S(n)*S(-n)
			 */
			
			IAtomList atomList = box.getLeafList();
			int numAtom = atomList.getAtomCount();
			double pl = 2*Math.PI/L;
			
			// S(n) and S(-n)
			
			double realMagnitude = 0.0;
			double imagMagnitude = 0.0;
			
			for (int i=0; i<numAtom; i++){
				IAtom atomi = atomList.getAtom(i);
				IVectorMutable posAtomi = ((IAtomPositioned)atomi).getPosition();
				double chargei = ((MyCharge)atomAgentManager.getAgent(atomi)).charge;
				
				double Sn = pl*(nVector[vecCounter].dot(posAtomi));
				
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
		IMoleculeList moleculeList = box.getMoleculeList();
		
		/*
		 * molecules can be monoatomic or multiatomic
		 */
		int numMolecules = moleculeList.getMoleculeCount();
		
		double prefcancel = alpha/Math.sqrt(Math.PI);
		double cancelTerm = 0.0;
		double pairTerm = 0.0;
		
		
		for (int i=0; i<numMolecules; i++){
			IMolecule molecule = moleculeList.getMolecule(i);
			
			/*
			 * Multi atomic
			 */
			
			IVectorMutable dabVector = space.makeVector();
			int numSites = molecule.getChildList().getAtomCount();
		
			// cancel-Term
			for (int a=0; a<numSites; a++){
				IAtom sitea = molecule.getChildList().getAtom(a);
				
				double chargea = ((MyCharge)atomAgentManager.getAgent(sitea)).charge;
				cancelTerm += prefcancel*chargea*chargea;
			}
			
			// self-Term
			moleculeBasis.atom = molecule;
			iterator.setBasis(moleculeBasis);
			iterator.reset();
			
			for (IAtomList pair = iterator.next(); pair!= null; pair = iterator.next()){
				IAtom sitea = pair.getAtom(0);
				IAtom siteb = pair.getAtom(1);
				
				IVectorMutable posAtoma = ((IAtomPositioned)sitea).getPosition();
				IVectorMutable posAtomb = ((IAtomPositioned)siteb).getPosition();
				
				double chargea = ((MyCharge)atomAgentManager.getAgent(sitea)).charge;
				double chargeb = ((MyCharge)atomAgentManager.getAgent(siteb)).charge;
			
				dabVector.Ev1Mv2(posAtomb, posAtoma);
				double dab = Math.sqrt(dabVector.squared());
				
				pairTerm += chargea*chargeb*(1-SpecialFunctions.erfc(alpha*dab))/dab;
			}
		
		}
		
		
		// taking the 0.5 out because the AtomPairIterator does not do the double counting
		return cancelTerm + pairTerm;
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
	
	protected final AtomLeafAgentManager atomAgentManager;
	protected final IVectorMutable[] nVector;
	protected final int nVectorMax;
	protected final double alpha;
	protected final IBox box;
	protected NeighborCriterion criterion;
	protected final AtomPair atomPair;
	protected AtomsetIteratorBasisDependent iterator;
	protected final MoleculeSetSinglet moleculeBasis;
	private final Space space;
	
	
	public static class MyCharge{
		
		public MyCharge(double charge){
			this.charge = charge;
		}
		
		public final double charge;
		
	}
}
