/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.species.ISpecies;

import java.util.LinkedList;

public class SimulationEventManager {


    private transient final LinkedList<SimulationListener> intervalListeners = new LinkedList<SimulationListener>();
    private final Simulation simulation;

    public SimulationEventManager(Simulation sim) {
        simulation = sim;
    }

    /* (non-Javadoc)
     * @see etomica.util.IEventManager#addListener(java.lang.Object)
     */
    public synchronized void addListener(SimulationListener listener) {
        if (listener == null) throw new NullPointerException("Cannot add null as a listener to Box");
//        if (listeners.contains(listener)) {
//            throw new RuntimeException(listener+" is already an interval action");
//        }
        intervalListeners.add(listener);
    }

    public synchronized void boxAdded(Box box) {
        SimulationBoxEvent e = new SimulationBoxEvent(simulation, box);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationBoxAdded(e);
        }
    }

    public synchronized void boxRemoved(Box box) {
        SimulationBoxEvent e = new SimulationBoxEvent(simulation, box);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationBoxRemoved(e);
        }
    }

    public synchronized void speciesAdded(ISpecies species) {
        SimulationSpeciesEvent e = new SimulationSpeciesEvent(simulation, species);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesAdded(e);
        }
    }

    public synchronized void speciesRemoved(ISpecies species) {
        SimulationSpeciesEvent e = new SimulationSpeciesEvent(simulation, species);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesRemoved(e);
        }
    }

    public synchronized void speciesIndexChanged(ISpecies species, int index) {
        SimulationSpeciesIndexEvent e = new SimulationSpeciesIndexEvent(simulation, species, index);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesIndexChanged(e);
        }
    }

    public synchronized void speciesMaxIndexChanged(int index) {
        SimulationIndexEvent e = new SimulationIndexEvent(simulation, index);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationSpeciesMaxIndexChanged(e);
        }
    }

    public synchronized void atomTypeIndexChanged(AtomType atomType, int index) {
        SimulationAtomTypeIndexEvent e = new SimulationAtomTypeIndexEvent(simulation, atomType, index);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationAtomTypeIndexChanged(e);
        }
    }

    public synchronized void atomTypeMaxIndexChanged(int index) {
        SimulationIndexEvent e = new SimulationIndexEvent(simulation, index);
        for (int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).simulationAtomTypeMaxIndexChanged(e);
        }
    }

    /* (non-Javadoc)
     * @see etomica.util.IEventManager#removeListener(java.lang.Object)
     */
    public synchronized void removeListener(SimulationListener listener) {
        intervalListeners.remove(listener);
    }
}
