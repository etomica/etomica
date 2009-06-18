package etomica.integrator;

import java.util.LinkedList;

import etomica.api.IIntegratorEventManager;
import etomica.api.IIntegratorListener;

public class IntegratorEventManager implements java.io.Serializable, IIntegratorEventManager {
    
    private final LinkedList<IIntegratorListener> intervalListeners = new LinkedList<IIntegratorListener>();


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
    
    public synchronized void stepStarted() {
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).integratorStepStarted(null);
        }
    }
    
    public synchronized void stepFinished() {
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).integratorStepFinished(null);
        }
        
    }
    
    public synchronized void initialized() {
        for(int i = 0; i < intervalListeners.size(); i++) {
            intervalListeners.get(i).integratorInitialized(null);
        }
        
    }
}

