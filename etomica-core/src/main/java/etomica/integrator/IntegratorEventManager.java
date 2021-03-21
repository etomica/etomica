/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import java.util.ArrayList;

public class IntegratorEventManager {

    protected final ArrayList<IntegratorListener> listeners = new ArrayList<IntegratorListener>();
    private final IntegratorEvent event;
    private boolean eventing;

    public IntegratorEventManager(Integrator integrator) {
        this.event = new IntegratorEvent(integrator);
    }

    public IntegratorListener[] getListeners() {
        return this.listeners.toArray(new IntegratorListener[]{});
    }

    /**
     * Adds the given listener to this event manager.
     */
    public void addListener(IntegratorListener newListener) {
        if (newListener == null) throw new NullPointerException("Cannot add null as a listener to Integrator");
        if (listeners.contains(newListener)) {
            throw new RuntimeException(newListener + " is already an interval action");
        }
        listeners.add(newListener);
    }

    /**
     * Removes the given listener from this event manager.
     */
    public void removeListener(IntegratorListener listener) {
        listeners.remove(listener);
    }

    /**
     * Returns true if the event manager is currently firing events.
     */
    public boolean firingEvent() {
        return eventing;
    }

    public void stepStarted() {
        eventing = true;
        for (IntegratorListener listener : listeners) {
            listener.integratorStepStarted(event);
        }
        eventing = false;
    }

    public void stepFinished() {
        eventing = true;
        for (IntegratorListener listener : listeners) {
            listener.integratorStepFinished(event);
        }
        eventing = false;
    }

    public void initialized() {
        eventing = true;
        for (IntegratorListener listener : listeners) {
            listener.integratorInitialized(event);
        }
        eventing = false;
    }

    public void forcePrecomputed() {
        eventing = true;
        for (IntegratorListener l : listeners) {
            if (l instanceof IntegratorListenerMD) {
                ((IntegratorListenerMD) l).integratorForcePrecomputed(event);
            }
        }
        eventing = false;
    }

    public void forceComputed() {
        eventing = true;
        for (IntegratorListener l : listeners) {
            if (l instanceof IntegratorListenerMD) {
                ((IntegratorListenerMD) l).integratorForceComputed(event);
            }
        }
        eventing = false;
    }

    public void preThermostat() {
        eventing = true;
        for (IntegratorListener l : listeners) {
            if (l instanceof IntegratorListenerMD) {
                ((IntegratorListenerMD) l).preThermostat(event);
            }
        }
        eventing = false;
    }

}

