package etomica.spin;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;


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
public class MCMoveSpinFlip extends MCMovePhase {

    /**
     * @param potentialMaster
     * @param nPhases
     */
    public MCMoveSpinFlip(PotentialMaster potentialMaster) {
        super(potentialMaster);
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
        AtomArrayList leafList = phase.getSpeciesMaster().leafList;
        atom = (AtomLeaf)leafList.get(Simulation.random.nextInt(leafList.size()));
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        atom.coord.position().TE(-1);
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
        atom.coord.position().TE(-1);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#affectedAtoms(etomica.Phase)
     */
    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#energyChange(etomica.Phase)
     */
    public double energyChange() {
        return uNew - uOld;
    }

    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected AtomLeaf atom;
    protected double uOld;
    protected double uNew = Double.NaN;

}
