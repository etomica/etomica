/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Convenience class used to define a SimulationAction. Implements all methods
 * of SimulationAction interface, except for actionPerformed.
 */
public abstract class SimulationActionAdapter implements SimulationAction, java.io.Serializable {

	/**
	 * @return Returns the simulation on which this action will be performed.
	 */
	public Simulation getSimulation() {
		return simulation;
	}

	/**
	 * @param simulation
	 *            The simulation on which this action will be performed.
	 */
	protected void setSimulation(Simulation simulation, Space _space) {
		this.simulation = simulation;
		this.space = _space;
	}

	protected Simulation simulation;
	protected Space space;
}
