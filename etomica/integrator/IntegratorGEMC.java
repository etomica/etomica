package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Space;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMoleculeExchange;
import etomica.integrator.mcmove.MCMoveVolumeExchange;
import etomica.potential.PotentialMaster;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator. Used to evaluate fluid-fluid
 * phase coexistence. Written to apply to only two phases.
 * 
 * @author David Kofke
 */

public class IntegratorGEMC extends IntegratorMC implements EtomicaElement {

    public Phase secondPhase;

    public IntegratorGEMC(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster);
        phaseCountMax = 2;
        phase = new Phase[phaseCountMax];
        atomDisplace0 = new MCMoveAtom(potentialMaster);
        atomDisplace1 = new MCMoveAtom(potentialMaster);
        volumeExchange = new MCMoveVolumeExchange(potentialMaster, space);
        moleculeExchange = new MCMoveMoleculeExchange(potentialMaster, space);
        setMCMoves(new MCMove[] { atomDisplace0, atomDisplace1, volumeExchange,
                moleculeExchange });
        atomDisplace0.setAdjustInterval(100);
        atomDisplace1.setAdjustInterval(100);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo(
                "Gibbs-ensemble Monte Carlo simulation of phase coexistence");
        return info;
    }

    public boolean addPhase(Phase p) {
        if(phase[0] == null) setPhase0(p);
        else if(phase[1] == null) setPhase1(p);
        else return false;
        return true;
    }

    private void setPhase0(Phase phase0) {
        if(phase0 == null) return;
        this.phase[0] = phase0;
        atomDisplace0.setPhase(new Phase[] {phase0});
        if(phase[1] != null) {
            volumeExchange.setPhase(phase);
            moleculeExchange.setPhase(phase);
        }
    }
    private void setPhase1(Phase phase1) {
        if(phase1 == null) return;
        this.phase[1] = phase1;
        atomDisplace1.setPhase(new Phase[] {phase1});
        if(phase[0] != null) {
            volumeExchange.setPhase(phase);
            moleculeExchange.setPhase(phase);
        }
    }

    public Phase[] getPhase() {
        return phase;
    }

    public void setPhase(Phase[] phases) {
        setPhase0(phases[0]);
        setPhase1(phases[1]);
    }

    /**
     * Returns the object that performs the volume-exchange move in the GE
     * simulation. Having handle to this object is needed to adjust trial
     * frequency and view acceptance rate.
     */
    public MCMoveVolumeExchange getMCMoveVolumeExchange() {
        return volumeExchange;
    }

    /**
     * Returns the object that performs the molecule-exchange move in the GE
     * simulation. Having handle to this object is needed to adjust trial
     * frequency and view acceptance rate.
     */
    public MCMoveMoleculeExchange getMCMoveMoleculeExchange() {
        return moleculeExchange;
    }

    /**
     * Returns the object that performs the atom-displacement move in the GE
     * simulation. Having handle to this object is needed to adjust trial
     * frequency and view acceptance rate.
     * 
     * @param i
     *            indicates request for move for 0th phase (i = 0) or first
     *            phase (i = 1)
     */
    public MCMoveAtom getMCMoveAtom(int i) {
        return (i == 0) ? atomDisplace0 : atomDisplace1;
    }

    private final MCMoveAtom atomDisplace0;
    private final MCMoveAtom atomDisplace1;
    private final MCMoveVolumeExchange volumeExchange;
    private final MCMoveMoleculeExchange moleculeExchange;

}