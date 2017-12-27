/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

public interface SimulationListener {

    default void simulationBoxAdded(SimulationBoxEvent e) {}

    default void simulationBoxRemoved(SimulationBoxEvent e) {}

    default void simulationSpeciesAdded(SimulationSpeciesEvent e) {}

    default void simulationSpeciesRemoved(SimulationSpeciesEvent e) {}

    default void simulationSpeciesIndexChanged(SimulationSpeciesIndexEvent e) {}

    default void simulationSpeciesMaxIndexChanged(SimulationIndexEvent e) {}

    default void simulationAtomTypeIndexChanged(SimulationAtomTypeIndexEvent e) {}

    default void simulationAtomTypeMaxIndexChanged(SimulationIndexEvent e) {}
}
