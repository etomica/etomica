package etomica.spin;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;


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

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#doTrial()
     */
    public boolean doTrial() {
        atom = phases[0].speciesMaster.atomList.getRandom();
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar(phases[0]);
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
        uNew = energyMeter.getDataAsScalar(phases[0]);
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
    public AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#energyChange(etomica.Phase)
     */
    public double energyChange(Phase phase) {
        return (this.phases[0] == phase) ? uNew - uOld : 0.0;
    }

    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected Atom atom;
    protected double uOld;
    protected double uNew = Double.NaN;

}
