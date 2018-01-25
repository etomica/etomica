/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.species.ISpecies;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.CopyOnWriteArrayList;

public class SimulationEventManager {

    private final List<SimulationListener> listeners = new ArrayList<>();
    private final Simulation simulation;

    protected SimulationEventManager(Simulation sim) {
        simulation = sim;
    }

    public void addListener(SimulationListener listener) {
        Objects.requireNonNull(listener);
        if (listeners.contains(listener)) {
            throw new IllegalArgumentException("This listener is already registered");
        }
        listeners.add(listener);
    }

    public void removeListener(SimulationListener listener) {
        listeners.remove(listener);
    }

    public void boxAdded(Box box) {
        SimulationBoxEvent e = new SimulationBoxEvent(simulation, box);
        for (SimulationListener listener : listeners) {
            listener.simulationBoxAdded(e);
        }
    }

    public void boxRemoved(Box box) {
        SimulationBoxEvent e = new SimulationBoxEvent(simulation, box);
        for (SimulationListener listener : listeners) {
            listener.simulationBoxRemoved(e);
        }
    }

    public void speciesAdded(ISpecies species) {
        SimulationSpeciesEvent e = new SimulationSpeciesEvent(simulation, species);
        for (SimulationListener listener : listeners) {
            listener.simulationSpeciesAdded(e);
        }
    }

    public void speciesRemoved(ISpecies species) {
        SimulationSpeciesEvent e = new SimulationSpeciesEvent(simulation, species);
        for (SimulationListener listener : listeners) {
            listener.simulationSpeciesRemoved(e);
        }
    }

    public void speciesIndexChanged(ISpecies species, int index) {
        SimulationSpeciesIndexEvent e = new SimulationSpeciesIndexEvent(simulation, species, index);
        for (SimulationListener listener : listeners) {
            listener.simulationSpeciesIndexChanged(e);
        }
    }

    public void speciesMaxIndexChanged(int index) {
        SimulationIndexEvent e = new SimulationIndexEvent(simulation, index);
        for (SimulationListener listener : listeners) {
            listener.simulationSpeciesMaxIndexChanged(e);
        }
    }

    public void atomTypeIndexChanged(AtomType atomType, int index) {
        SimulationAtomTypeIndexEvent e = new SimulationAtomTypeIndexEvent(simulation, atomType, index);
        for (SimulationListener listener : listeners) {
            listener.simulationAtomTypeIndexChanged(e);
        }
    }

    public void atomTypeMaxIndexChanged(int index) {
        SimulationIndexEvent e = new SimulationIndexEvent(simulation, index);
        for (SimulationListener listener : listeners) {
            listener.simulationAtomTypeMaxIndexChanged(e);
        }
    }
}
