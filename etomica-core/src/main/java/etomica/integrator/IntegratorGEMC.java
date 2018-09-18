/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.util.random.IRandom;
import etomica.integrator.mcmove.MCMoveMoleculeExchange;
import etomica.integrator.mcmove.MCMoveVolumeExchange;
import etomica.space.Space;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator. Used to evaluate fluid-fluid
 * box coexistence. Written to apply to only two boxs.
 * 
 * @author David Kofke
 */
public class IntegratorGEMC extends IntegratorManagerMC {

    public IntegratorGEMC(IRandom random, Space space) {
        super(random);
        this.space = space;
    }

    public void addIntegrator(Integrator newIntegrator) {
        if (!(newIntegrator instanceof IntegratorBox)) {
            throw new IllegalArgumentException("Sub integrators must be able to handle a box");
        }
        if (integrators.size() == 2) {
            throw new IllegalArgumentException("Only 2 sub-integrators can be added");
        }
        super.addIntegrator(newIntegrator);
        if (integrators.size() == 2) {
            volumeExchange = new MCMoveVolumeExchange(((IntegratorBox)newIntegrator).getPotentialMaster(), random,
                    space, (IntegratorBox)integrators.get(0),(IntegratorBox) integrators.get(1));
            moleculeExchange = new MCMoveMoleculeExchange(((IntegratorBox)newIntegrator).getPotentialMaster(), random,
                    space, (IntegratorBox) integrators.get(0),(IntegratorBox) integrators.get(1));
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
    private Space space;

}
