/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import java.util.ArrayList;

import etomica.api.IIntegratorEventManager;
import etomica.api.IIntegratorListener;

public class IntegratorEventManager implements IIntegratorEventManager {
    
    private final ArrayList<IIntegratorListener> intervalListeners = new ArrayList<IIntegratorListener>();
    protected boolean eventing;


    public synchronized void addListener(IIntegratorListener newListener) {
        if(newListener == null) throw new NullPointerException("Cannot add null as a listener to Integrator");
        if (intervalListeners.contains(newListener)) {
            throw new RuntimeException(newListener+" is already an interval action");
        }
        intervalListeners.add(newListener);
    }

    public synchronized void removeListener(IIntegratorListener listener) {
        intervalListeners.remove(listener);
    }
    
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
}

