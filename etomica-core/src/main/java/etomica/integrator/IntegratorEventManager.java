/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.api.IIntegratorListener;
import etomica.api.IIntegratorListenerMD;

import java.util.ArrayList;

public class IntegratorEventManager {

    protected final ArrayList<IIntegratorListener> intervalListeners = new ArrayList<IIntegratorListener>();
    protected boolean eventing;

    /**
     * Adds the given listener to this event manager.
     */
    public synchronized void addListener(IIntegratorListener newListener) {
        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Integrator");
        if (intervalListeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        intervalListeners.add(newListener);
    }

    /**
     * Removes the given listener from this event manager.
     */
    public synchronized void removeListener(IIntegratorListener listener) {
        intervalListeners.remove(listener);
    }

    /**
     * Returns true if the event manager is currently firing events.
     */
    public boolean firingEvent() {
        return eventing;
    }

    public synchronized void stepStarted() {
        eventing = true;
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).integratorStepStarted(null);
        }
        eventing = false;
    }

    public synchronized void stepFinished() {
        eventing = true;
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).integratorStepFinished(null);
        }
        eventing = false;
    }

    public synchronized void initialized() {
        eventing = true;
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).integratorInitialized(null);
        }
        eventing = false;
    }
    
    public synchronized void forcePrecomputed() {
        eventing = true;
        for(int i = 0; i < intervalListeners.size(); i++) {
            IIntegratorListener l = intervalListeners.get(i);
            if (l instanceof IIntegratorListenerMD) {
                ((IIntegratorListenerMD)intervalListeners.get(i)).integratorForcePrecomputed(null);
            }
        }
        eventing = false;
    }

    public synchronized void forceComputed() {
        eventing = true;
        for(int i = 0; i < intervalListeners.size(); i++) {
            IIntegratorListener l = intervalListeners.get(i);
            if (l instanceof IIntegratorListenerMD) {
                ((IIntegratorListenerMD)intervalListeners.get(i)).integratorForceComputed(null);
            }
        }
        eventing = false;
    }
}

