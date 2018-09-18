/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.box.Box;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

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
        for (SimulationListener listener : listeners) {
            listener.simulationBoxAdded(simulation, box);
        }
    }

    public void boxRemoved(Box box) {
        for (SimulationListener listener : listeners) {
            listener.simulationBoxRemoved(simulation, box);
        }
    }
}
