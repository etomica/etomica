/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.listener;

import etomica.integrator.IntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.util.Arrays;

public class IntegratorListenerGroupSeries implements IIntegratorListener, java.io.Serializable {

    /**
     * Constructs an action group that holds no actions.
     */
    public IntegratorListenerGroupSeries() {
        this(new IIntegratorListener[0]);
    }
    
    /**
     * Defines group via the given array of actions.  Copy
     * of array is made and used internally.
     */
    public IntegratorListenerGroupSeries(IIntegratorListener[] listeners) {
        this.listeners = listeners.clone();
        intervalCount = 0;
        interval = 1;
    }
    
    public void integratorInitialized(IntegratorEvent e) {
        for(int i=0; i<listeners.length; i++) {
            listeners[i].integratorInitialized(e);
        }
    }
    
    public void integratorStepStarted(IntegratorEvent e) {
        intervalCount++;
        if(intervalCount >= interval) {
            for(int i=0; i<listeners.length; i++) {
                listeners[i].integratorStepStarted(e);
            }
        }
    }
    
    public void integratorStepFinished(IntegratorEvent e) {
        if(intervalCount >= interval) {
            for(int i=0; i<listeners.length; i++) {
                listeners[i].integratorStepFinished(e);
            }
            intervalCount = 0;
        }
    }
    
    /**
     * Adds the given action to the group.  No check is made of whether
     * action is already in group; it is added regardless.  
     * @param newListener
     */
    public void addListener(IIntegratorListener newListener) {
        listeners = (IIntegratorListener[])Arrays.addObject(listeners, newListener);
    }
    
    /**
     * Removes the given action from the group.  No warning or
     * error is given if action is not in the group already.
     */
    public boolean removeListener(IIntegratorListener oldListener) {
        int num = listeners.length;
        listeners = (IIntegratorListener[])Arrays.removeObject(listeners, oldListener);
        return listeners.length != num; 
    }
    
    public IIntegratorListener[] getAllListeners() {
        return listeners.clone();
    }
    
    public void setInterval(int i) {
        interval = i;
    }
    
    private static final long serialVersionUID = 1L;
    private int interval;
    private int intervalCount;
    private IIntegratorListener[] listeners;
}
