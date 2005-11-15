package etomica.spin;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 22, 2005 by kofke
 */
public class MCMoveSpinFlip extends MCMove {

    /**
     * @param potentialMaster
     * @param nPhases
     */
    public MCMoveSpinFlip(PotentialMaster potentialMaster) {
        super(potentialMaster, 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        setName("MCMoveSpinFlip");

    }

    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
    }
    
    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#doTrial()
     */
    public boolean doTrial() {
        atom = phase.getSpeciesMaster().atomList.getRandom();
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        atom.coord.position().TE(-1);
        uNew = Double.NaN;
        return true;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#lnTrialRatio()
     */
    public double lnTrialRatio() {
        return 0.0;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#lnProbabilityRatio()
     */
    public double lnProbabilityRatio() {
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld)/temperature;
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
        atom.coord.position().TE(-1);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#affectedAtoms(etomica.Phase)
     */
    public AtomIterator affectedAtoms(Phase p) {
        if(p != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#energyChange(etomica.Phase)
     */
    public double energyChange(Phase p) {
        return (p == phase) ? uNew - uOld : 0.0;
    }

    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected Atom atom;
    protected double uOld;
    protected double uNew = Double.NaN;

}
