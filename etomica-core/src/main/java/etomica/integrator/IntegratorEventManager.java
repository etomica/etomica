/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.api.IIntegratorListener;
import etomica.api.IIntegratorListenerMD;

import java.util.ArrayList;

public class IntegratorEventManager {

    protected final ArrayList<IIntegratorListener> listeners = new ArrayList<IIntegratorListener>();
    private final IntegratorEvent event;
    protected boolean eventing;

    public IntegratorEventManager(Integrator integrator) {
        this.event = new IntegratorEvent(integrator);
    }

    /**
     * Adds the given listener to this event manager.
     */
    public synchronized void addListener(IIntegratorListener newListener) {
        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Integrator");
        if (listeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        listeners.add(newListener);
    }

    /**
     * Removes the given listener from this event manager.
     */
    public synchronized void removeListener(IIntegratorListener listener) {
        listeners.remove(listener);
    }

    /**
     * Returns true if the event manager is currently firing events.
     */
    public boolean firingEvent() {
        return eventing;
    }

    public synchronized void stepStarted() {
        eventing = true;
        for (IIntegratorListener listener : listeners) {
            listener.integratorStepStarted(event);
        }
        eventing = false;
    }

    public synchronized void stepFinished() {
        eventing = true;
        for (IIntegratorListener listener : listeners) {
            listener.integratorStepFinished(event);
        }
        eventing = false;
    }

    public synchronized void initialized() {
        eventing = true;
        for (IIntegratorListener listener : listeners) {
            listener.integratorInitialized(event);
        }
        eventing = false;
    }
    
    public synchronized void forcePrecomputed() {
        eventing = true;
        for (IIntegratorListener l : listeners) {
            if (l instanceof IIntegratorListenerMD) {
                ((IIntegratorListenerMD) l).integratorForcePrecomputed(event);
            }
        }
        eventing = false;
    }

    public synchronized void forceComputed() {
        eventing = true;
        for (IIntegratorListener l : listeners) {
            if (l instanceof IIntegratorListenerMD) {
                ((IIntegratorListenerMD) l).integratorForceComputed(event);
            }
        }
        eventing = false;
    }
}

