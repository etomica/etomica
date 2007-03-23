package etomica.integrator.mcmove;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomSource;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateAngular;
import etomica.space.Orientation;
import etomica.util.IRandom;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */
public class MCMoveRotate extends MCMovePhaseStep {
    
    private static final long serialVersionUID = 2L;
    private MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final Orientation oldOrientation;
    private AtomSource atomSource;

    private transient AtomLeaf molecule;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient Orientation orientation;
    
    protected final IRandom random;

    public MCMoveRotate(Simulation sim) {
        this(sim.getPotentialMaster(), sim.getRandom());
    }
    
    public MCMoveRotate(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        oldOrientation = potentialMaster.getSpace().makeOrientation();
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }
    
    /**
     * Sets the AtomSource used to select Atoms acted on by MC trials.
     */
    public void setAtomSource(AtomSource newAtomSource) {
        atomSource = newAtomSource;
    }
     
    /**
     * Returns the AtomSource used to select Atoms acted on by MC trials.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
     
    public boolean doTrial() {
        if(phase.moleculeCount()==0) {return false;}
        molecule = (AtomLeaf)atomSource.getAtom();

        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        orientation = ((ICoordinateAngular)molecule.getCoord()).getOrientation(); 
        oldOrientation.E(orientation);  //save old orientation
        orientation.randomRotation(random, stepSize);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {return 1.0;}
    
    public double getB() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        orientation.E(oldOrientation);
    }

    public double energyChange() {return uNew - uOld;}
    
    public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
    }
}
