/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.osmoticvirial;

import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator. Used to evaluate fluid-fluid
 * box coexistence. Written to apply to only two boxs.
 * 
 * @author David Kofke
 */
public class IntegratorRGEMC extends IntegratorManagerMC {

    private MCMoveGeometricClusterRestrictedGE mcMoveGeometricClusterRestrictedGE;
    private Space space;
    private final ISpecies seed;

    public IntegratorRGEMC(IRandom random, Space space, ISpecies seed) {
        super(random);
        this.space = space;
        this.seed = seed;
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

            mcMoveGeometricClusterRestrictedGE =
                    new MCMoveGeometricClusterRestrictedGE(((IntegratorBox) newIntegrator).getPotentialMaster(),
                            space, random, ((IntegratorBox) integrators.get(0)).getBox(), ((IntegratorBox) integrators.get(1)).getBox(), seed);
            moveManager.recomputeMoveFrequencies();
            moveManager.addMCMove(mcMoveGeometricClusterRestrictedGE);
        }
    }

    /**
     * Returns the object that performs the molecule-exchange move in the GE
     * simulation. Having handle to this object is needed to adjust trial
     * frequency and view acceptance rate.
     */
    public MCMoveGeometricClusterRestrictedGE getMCMoveMoleculeExchange() {
        return mcMoveGeometricClusterRestrictedGE;
    }

}
