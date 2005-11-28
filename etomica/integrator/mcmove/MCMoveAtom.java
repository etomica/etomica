package etomica.integrator.mcmove;

import etomica.atom.Atom;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/09/02 Added energyChange() method 
  * 07/09/03 Changed fields from private to protected
  */
  
public class MCMoveAtom extends MCMove {
    
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected Atom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;

    public MCMoveAtom(Simulation sim) {
        this(sim.potentialMaster, sim.getDefaults().atomSize, sim.getDefaults().boxSize/2, 
                sim.getDefaults().ignoreOverlap);
    }
    
    public MCMoveAtom(PotentialMaster potentialMaster, double stepSize, double stepSizeMax,
            boolean fixOverlap) {
        super(potentialMaster, 1);
        atomSource = new AtomSourceRandomLeaf();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = potentialMaster.getSpace().makeVector();
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        setName("MCMoveAtom");
        this.fixOverlap = fixOverlap;
    }
    
    public final Dimension getStepSizeDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMaxDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMinDimension() {return Dimension.LENGTH;}
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        if(uOld > 1e10 && !fixOverlap) {
            throw new RuntimeException(new ConfigurationOverlapException(atom.node.parentPhase()));
        }
        translationVector.setRandomCube();
        translationVector.TE(stepSize);
        atom.coord.position().PE(translationVector);
        uNew = Double.NaN;
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
//        energyMeter.setTarget(atom);
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public double energyChange(Phase p) {return uNew - uOld;}
    
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
        translationVector.TE(-1);
        atom.coord.position().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms(Phase p) {
        if(phase != p) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
        atomSource.setPhase(p);
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