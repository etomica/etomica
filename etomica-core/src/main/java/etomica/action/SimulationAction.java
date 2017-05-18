/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.simulation.Simulation;

/**
 * Interface for classes that apply some elementary action (transformation) to a
 * simulation.
 *  
 */
public interface SimulationAction extends IAction {

	public Simulation getSimulation();
}
