package etomica.normalmode;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.box.Box;
import etomica.potential.Potential2;
import etomica.potential.PotentialMaster;
import etomica.space.IVectorRandom;
import etomica.util.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move.  Two atoms are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomCoupled extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final MeterPotentialEnergy energyMeter;
    protected final IVectorRandom translationVector;
    protected IAtom atom0, atom1;
    protected double uOld;
    protected double uNew;
    protected AtomSource atomSource;
    protected final IRandom random;
    protected Potential2 potential;
    protected final AtomPair pair;

    public MCMoveAtomCoupled(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = (IVectorRandom)potentialMaster.getSpace().makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        pair = new AtomPair();
    }
    
    public void setPotential(Potential2 newPotential) {
        potential = newPotential;
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom0 = atomSource.getAtom();
        atom1 = atomSource.getAtom();
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
        ((IAtomPositioned)atom0).getPosition().PE(translationVector);
        ((IAtomPositioned)atom1).getPosition().ME(translationVector);

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
        ((IAtomPositioned)atom0).getPosition().ME(translationVector);
        ((IAtomPositioned)atom1).getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom0);
        affectedAtomList.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
        potential.setBox(p);
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }
}