/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.util.ArrayList;
import java.util.List;

import etomica.action.IAction;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorHard;
import etomica.space.Space;

/**
 * Custom DMD integrator that handles hybrid simulations with
 * insertion/deletions.
 *  
 * @author Andrew Schultz
 */
public class IntegratorHardMDMC extends IntegratorHard {
    protected List<IAction> thermostatActions;

    public IntegratorHardMDMC(Simulation sim, PotentialMaster potentialMaster, Space _space, Box box) {
        super(sim, potentialMaster, _space, box);
        thermostatActions = new ArrayList<IAction>();
    }
    
    public void setBox(Box box) {
        super.setBox(box);
    }
    
    public void addThermostatAction(IAction a) {
        thermostatActions.add(a);
    }

    public void reset() {
        if (!thermostatting) {
            super.reset();
            return;
        }

        // we get here because HYBRID_MC thermostat calls reset.
        // neighbors should be fine, we just need to update the potential energy
        // IntegratorHard would reset collision times.  that will happen anyway from randomizeMomenta
        currentPotentialEnergy = meterPE.getDataAsScalar();
        if (currentPotentialEnergy == Double.POSITIVE_INFINITY) {
            System.err.println("overlap in configuration for "+box+" when resetting integrator");
            throw new ConfigurationOverlapException(box);
        }
        currentKineticEnergy = meterKE.getDataAsScalar();

    }
    
    public void doThermostat() {
        thermostatting = true;
        for (IAction a : thermostatActions) {
            a.actionPerformed();
        }
        super.doThermostat();
        thermostatting = false;
    }
}
