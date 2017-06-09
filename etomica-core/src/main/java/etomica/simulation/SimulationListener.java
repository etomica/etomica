/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

public interface SimulationListener {

    void simulationBoxAdded(SimulationBoxEvent e);

    void simulationBoxRemoved(SimulationBoxEvent e);

    void simulationSpeciesAdded(SimulationSpeciesEvent e);

    void simulationSpeciesRemoved(SimulationSpeciesEvent e);

    void simulationSpeciesIndexChanged(SimulationSpeciesIndexEvent e);

    void simulationSpeciesMaxIndexChanged(SimulationIndexEvent e);

    void simulationAtomTypeIndexChanged(SimulationAtomTypeEvent e);

    void simulationAtomTypeMaxIndexChanged(SimulationIndexEvent e);
}
