package etomica.spin;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.util.IRandom;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class MCMoveSpinFlip extends MCMoveBox {

    /**
     * @param potentialMaster
     * @param nBoxs
     */
    public MCMoveSpinFlip(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
    }
    
    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#doTrial()
     */
    public boolean doTrial() {
        AtomSet leafList = box.getSpeciesMaster().getLeafList();
        atom = (IAtomPositioned)leafList.getAtom(random.nextInt(leafList.getAtomCount()));
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        atom.getPosition().TE(-1);
        uNew = Double.NaN;
        return true;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#lnTrialRatio()
     */
    public double getA() {
        return 1.0;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#lnProbabilityRatio()
     */
    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#acceptNotify()
     */
    public void acceptNotify() {
        //nothing to do
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#rejectNotify()
     */
    public void rejectNotify() {
        atom.getPosition().TE(-1);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#affectedAtoms(etomica.Box)
     */
    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#energyChange(etomica.Box)
     */
    public double energyChange() {
        return uNew - uOld;
    }

    private static final long serialVersionUID = 1L;
    protected final IRandom random;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected IAtomPositioned atom;
    protected double uOld;
    protected double uNew = Double.NaN;
}
