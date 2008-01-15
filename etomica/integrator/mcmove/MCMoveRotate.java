package etomica.integrator.mcmove;

import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.IAtomOriented;
import etomica.atom.IMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialMaster;
import etomica.space.IOrientation;
import etomica.util.IRandom;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */
public class MCMoveRotate extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    private MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final IOrientation oldOrientation;
    private AtomSource atomSource;

    private transient IAtomOriented molecule;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient IOrientation iOrientation;
    
    protected final IRandom random;

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
        atomSource = new AtomSourceRandomMolecule();
        ((AtomSourceRandomMolecule)atomSource).setRandom(random);
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
        if(box.moleculeCount()==0) {return false;}
        molecule = (IAtomOriented)((IMolecule)atomSource.getAtom()).getChildList().getAtom(0);

        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        iOrientation = molecule.getOrientation(); 
        oldOrientation.E(iOrientation);  //save old orientation
        iOrientation.randomRotation(random, stepSize);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }//end of doTrial
    
    public double getA() {return 1.0;}
    
    public double getB() {
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        iOrientation.E(oldOrientation);
    }

    public double energyChange() {return uNew - uOld;}
    
    public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
    }
}
