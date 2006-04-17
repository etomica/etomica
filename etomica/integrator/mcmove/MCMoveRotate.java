package etomica.integrator.mcmove;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.ICoordinateAngular;
import etomica.space.Orientation;
import etomica.space.Space;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */

 /* History of changes
  * 7/9/02 Added energyChange() method
  */
  
public class MCMoveRotate extends MCMovePhaseStep {
    
    private MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final Orientation oldOrientation;

    private transient AtomLeaf molecule;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient Orientation orientation;

    public MCMoveRotate(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        oldOrientation = space.makeOrientation();
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }
     
    public boolean doTrial() {
        if(phase.moleculeCount()==0) {return false;}
        molecule = (AtomLeaf)phase.randomMolecule();

        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        orientation = ((ICoordinateAngular)molecule.coord).orientation(); 
        oldOrientation.E(orientation);  //save old orientation
        orientation.randomRotation(stepSize);
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
