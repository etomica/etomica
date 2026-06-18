/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.osmoticvirial;

import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.potential.IPotential2;
import etomica.potential.compute.NeighborManager;
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
    private final NeighborManager[] neighborManagers = new NeighborManager[2];
    private final IPotential2[][] pairPotentials;

    public IntegratorRGEMC(IRandom random, Space space, ISpecies seed, IPotential2[][] pairPotentials) {
        super(random);
        this.space = space;
        this.seed = seed;
        this.pairPotentials = pairPotentials;
    }

    public void addIntegrator(IntegratorBox newIntegrator, NeighborManager nbrManager) {
        neighborManagers[integrators.size()] = nbrManager;
        if (integrators.size() == 2) {
            throw new IllegalArgumentException("Only 2 sub-integrators can be added");
        }
        super.addIntegrator(newIntegrator);
        if (integrators.size() == 2) {

            mcMoveGeometricClusterRestrictedGE =
                    new MCMoveGeometricClusterRestrictedGE(neighborManagers[0], neighborManagers[1],
                            space, random, ((IntegratorBox) integrators.get(0)).getBox(), ((IntegratorBox) integrators.get(1)).getBox(), seed, pairPotentials);
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
