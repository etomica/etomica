/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.integrator.mcmove.MCMoveMoleculeExchangeFasterer;
import etomica.integrator.mcmove.MCMoveVolumeExchangeFasterer;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator. Used to evaluate fluid-fluid
 * box coexistence. Written to apply to only two boxs.
 *
 * @author David Kofke
 */
public class IntegratorGEMC {

    public static IntegratorManagerMC buildGEMC(IntegratorBoxFasterer integrator1, IntegratorBoxFasterer integrator2, IRandom random, Space space) {
        IntegratorManagerMC integratorGEMC = new IntegratorManagerMC(random);
        integratorGEMC.addIntegrator(integrator1);
        integratorGEMC.addIntegrator(integrator2);
        MCMoveVolumeExchangeFasterer volumeExchange = new MCMoveVolumeExchangeFasterer(random, space, integrator1, integrator2);
        MCMoveMoleculeExchangeFasterer moleculeExchange = new MCMoveMoleculeExchangeFasterer(random, space, integrator1, integrator2);
        integratorGEMC.getMoveManager().recomputeMoveFrequencies();
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        return integratorGEMC;
    }
}
