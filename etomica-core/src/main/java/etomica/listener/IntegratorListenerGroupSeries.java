/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.listener;

import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.util.Arrays;

public class IntegratorListenerGroupSeries implements IntegratorListener, java.io.Serializable {

    /**
     * Constructs an action group that holds no actions.
     */
    public IntegratorListenerGroupSeries() {
        this(new IntegratorListener[0]);
    }
    
    /**
     * Defines group via the given array of actions.  Copy
     * of array is made and used internally.
     */
    public IntegratorListenerGroupSeries(IntegratorListener[] listeners) {
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
    public void addListener(IntegratorListener newListener) {
        listeners = (IntegratorListener[])Arrays.addObject(listeners, newListener);
    }
    
    /**
     * Removes the given action from the group.  No warning or
     * error is given if action is not in the group already.
     */
    public boolean removeListener(IntegratorListener oldListener) {
        int num = listeners.length;
        listeners = (IntegratorListener[])Arrays.removeObject(listeners, oldListener);
        return listeners.length != num; 
    }
    
    public IntegratorListener[] getAllListeners() {
        return listeners.clone();
    }
    
    public void setInterval(int i) {
        interval = i;
    }
    
    private static final long serialVersionUID = 1L;
    private int interval;
    private int intervalCount;
    private IntegratorListener[] listeners;
}
