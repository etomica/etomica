/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.integrator.mcmove.MCMoveMoleculeExchange;
import etomica.integrator.mcmove.MCMoveVolumeExchange;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator. Used to evaluate fluid-fluid
 * box coexistence. Written to apply to only two boxs.
 *
 * @author David Kofke
 */
public class IntegratorGEMC {

    public static IntegratorManagerMC buildGEMC(IntegratorBox integrator1, IntegratorBox integrator2, IRandom random, Space space) {
        IntegratorManagerMC integratorGEMC = new IntegratorManagerMC(random);
        integratorGEMC.setTemperature(integrator1.getTemperature());
        integratorGEMC.addIntegrator(integrator1);
        integratorGEMC.addIntegrator(integrator2);
        MCMoveVolumeExchange volumeExchange = new MCMoveVolumeExchange(random, space, integrator1, integrator2);
        MCMoveMoleculeExchange moleculeExchange = new MCMoveMoleculeExchange(random, space, integrator1, integrator2);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 1);
        int n1 = integrator1.getBox().getMoleculeList().size();
        int n2 = integrator2.getBox().getMoleculeList().size();
        integratorGEMC.getMoveManager().setFrequency(moleculeExchange, n1+n2);
        return integratorGEMC;
    }
}
