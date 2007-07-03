package etomica.integrator;

import etomica.EtomicaInfo;
import etomica.integrator.mcmove.MCMoveMoleculeExchange;
import etomica.integrator.mcmove.MCMoveVolumeExchange;
import etomica.util.IRandom;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator. Used to evaluate fluid-fluid
 * box coexistence. Written to apply to only two boxs.
 * 
 * @author David Kofke
 */
public class IntegratorGEMC extends IntegratorManagerMC {

    public IntegratorGEMC(IRandom random) {
        super(random);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo(
                "Gibbs-ensemble Monte Carlo simulation of box coexistence");
        return info;
    }

    public void addIntegrator(IIntegrator newIntegrator) {
        if (!(newIntegrator instanceof IntegratorBox)) {
            throw new IllegalArgumentException("Sub integrators must be able to handle a box");
        }
        if (nIntegrators == 2) {
            throw new IllegalArgumentException("Only 2 sub-integrators can be added");
        }
        super.addIntegrator(newIntegrator);
        if (nIntegrators == 2) {
            volumeExchange = new MCMoveVolumeExchange(((IntegratorBox)newIntegrator).getPotential(), random,
                    (IntegratorBox)integrators[0],(IntegratorBox)integrators[1]);
            moleculeExchange = new MCMoveMoleculeExchange(((IntegratorBox)newIntegrator).getPotential(), random,
                    (IntegratorBox)integrators[0],(IntegratorBox)integrators[1]);
            moveManager.recomputeMoveFrequencies();
            moveManager.addMCMove(volumeExchange);
            moveManager.addMCMove(moleculeExchange);
        }
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

    private static final long serialVersionUID = 1L;
    private MCMoveVolumeExchange volumeExchange;
    private MCMoveMoleculeExchange moleculeExchange;

}