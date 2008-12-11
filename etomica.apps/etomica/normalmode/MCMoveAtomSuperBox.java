 package etomica.normalmode;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.Potential2;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Standard Monte Carlo atom-displacement trial move. Two atoms are moved at a
 * time in such a way that all the atoms' "images" are being moved at the same time too.
 * 
 * This class is written specifically for 32 system size FCC packing
 * 
 * 
 * @author Tai Boon Tan
 */
public class MCMoveAtomSuperBox extends MCMoveBoxStep {


    public MCMoveAtomSuperBox(IPotentialMaster potentialMaster, IRandom random,
    		                 ISpace _space, CoordinateDefinition coordinateDefinition) {
        super(potentialMaster);
        this.random = random;
        this.coordinateDefinition = coordinateDefinition;
 
        cells = coordinateDefinition.getBasisCells();
        
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = (IVectorRandom)_space.makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        pair = new AtomPair();
        
      
        /*
         * determine the 'imaginary' box index
         * 	these boxes are not the actual simulation box!!!
         * 
         * boxIndex[a] = x
         * 	a is the cell # [total number of cells is 3x3x3x8 =216]
         * 	x is the box # [total number of boxes is 3x3x3 = 27]
         */
        boxIndex = new int[216];
        boxCells = new BasisCell[27][8];
    	int[] n = new int[27];
    	for (int j=0; j<6; j++){
    	
        	for(int i=0; i<cells.length/6; i++){
        		
        		int thisBoxIndex =  (i/12)*3 + ((i%6)/2) + ((j/2)*9);
        		boxIndex[i+(j*36)] = thisBoxIndex;
        		boxCells[thisBoxIndex][n[thisBoxIndex]] = cells[i+(j*36)];
        		n[thisBoxIndex]++;
        	}
        	
    	}
    	
        
    }
    
    public void setPotential(Potential2 newPotential) {
        potential = newPotential;
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {

    	randomNumber0 = random.nextInt(8);
    	BasisCell randomCell0 = boxCells[0][randomNumber0];
    	randomAtom0 = random.nextInt(4);
    	atom0 = randomCell0.molecules.getMolecule(randomAtom0).getChildList().getAtom(0);
    	
    	randomNumber1 = random.nextInt(8);
    	BasisCell randomCell1 = boxCells[0][randomNumber1];
    	randomAtom1 = random.nextInt(4);
    	atom1 = randomCell1.molecules.getMolecule(randomAtom1).getChildList().getAtom(0);
        
        if (atom0 == null || atom1 == null || atom0 == atom1) return false;
        energyMeter.setTarget(atom0);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom1);
        uOld += energyMeter.getDataAsScalar();
        pair.atom0 = atom0;
        pair.atom1 = atom1;
        
        uOld -= potential.energy(pair);
        if(uOld > 1e10) {
            throw new RuntimeException(new ConfigurationOverlapException(box));
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        
        /*
         * Moving all the corresponding atoms in other boxes
         * 	at the same time with the atom that was randomly 
         * 	generated 
         * 
         * there are total of 27 boxes
         */
        
        
        for (int i=0; i<27;i++){
        
        	((IAtomPositioned)boxCells[i][randomNumber0].molecules.getMolecule(randomAtom0).getChildList().getAtom(0)).getPosition().PE(translationVector);
        	((IAtomPositioned)boxCells[i][randomNumber1].molecules.getMolecule(randomAtom1).getChildList().getAtom(0)).getPosition().ME(translationVector);
 
        }
        
        uNew = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom0);
        uNew += energyMeter.getDataAsScalar();
        uNew -= potential.energy(pair);
        return true;
    }//end of doTrial
    
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {return 1.0;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        return -(uNew - uOld);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
    	
    	//System.err.println("Rejected atom is in cell["+randomNumber0+"]["+randomAtom0+"]");
    	
        for (int i=0; i<27;i++){
            
        	((IAtomPositioned)boxCells[i][randomNumber0].molecules.getMolecule(randomAtom0).getChildList().getAtom(0)).getPosition().ME(translationVector);
        	((IAtomPositioned)boxCells[i][randomNumber1].molecules.getMolecule(randomAtom1).getChildList().getAtom(0)).getPosition().PE(translationVector);
 
        }
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom0);
        affectedAtomList.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        potential.setBox(p);
    }
  
    
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final MeterPotentialEnergy energyMeter;
    protected final IVectorRandom translationVector;
    protected IAtomLeaf atom0, atom1;
    protected double uOld;
    protected double uNew;
    protected final IRandom random;
    protected Potential2 potential;
    protected final AtomPair pair;
    protected CoordinateDefinition coordinateDefinition;
    protected BasisCell[] cells;
    protected int[] boxIndex;
    protected BasisCell[][] boxCells; 
    protected int randomNumber0, randomNumber1, randomAtom0, randomAtom1;
}