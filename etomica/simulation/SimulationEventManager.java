package etomica.simulation;
import java.io.Serializable;

import etomica.action.activity.ControllerEvent;
import etomica.action.activity.ControllerListener;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveListener;
import etomica.lattice.LatticeEvent;
import etomica.lattice.LatticeListener;
import etomica.phase.PhaseEvent;
import etomica.phase.PhaseListener;

/**
 * Class to take care of listener lists and event firing for simulation elements.
 * A class can make an instance of this manager as a field, and delegate any 
 * listener management functions to it.  
 *
 * @see SimulationEvent
 * @see SimulationListener
 */
public class SimulationEventManager implements java.io.Serializable {
    
    /**
     * Adds a listener.  Synchronized to avoid conflict with removeListener.
     */
    public synchronized void addListener(Object listener) {
        //add listener to beginning of list 
        //placement at end causes problem if a listener removes and then adds itself to the list as part of its response to the event
        first = new SimulationEventManager.Linker(listener, first, null);
        if(first.next != null) first.next.previous = first;
        listenerCount++;
    }

    /**
     * Removes a listener.  Synchronized to avoid conflict with addListener.
     */
    public synchronized void removeListener(Object listener) {
        for(SimulationEventManager.Linker link=first; link!=null; link=link.next) {
            if(link.listener == listener) {
            	listenerCount--;
                if(link == first) {first = link.next;}
                if(link.previous != null) link.previous.next = link.next;
                if(link.next != null) link.next.previous = link.previous;
                return;
            }
        }
    }
    
    /**
     * @return the number of listeners currently managed by this instance.
     */
    public int listenerCount() {
    	return listenerCount;
    }
    
    /**
     * Returns the first of the linked-list of listeners.  Useful if 
     * using this class just to maintain list of listeners, while firing
     * events to them using methods outside this class.
     */
     public SimulationEventManager.Linker first() {return first;}

    public void fireEvent(PhaseEvent event) {
        for(SimulationEventManager.Linker link=first; link!=null; link=link.next) {
            ((PhaseListener)link.listener).actionPerformed(event);
        }
    }

    public void fireEvent(ControllerEvent event) {
        for(SimulationEventManager.Linker link=first; link!=null; link=link.next) {
            ((ControllerListener)link.listener).actionPerformed(event);
        }
    }

    public void fireEvent(MCMoveEvent event) {
        for(SimulationEventManager.Linker link=first; link!=null; link=link.next) {
            ((MCMoveListener)link.listener).actionPerformed(event);
        }
    }

    public void fireEvent(LatticeEvent event) {
        for(SimulationEventManager.Linker link=first; link!=null; link=link.next) {
            ((LatticeListener)link.listener).actionPerformed(event);
        }
    }

    private SimulationEventManager.Linker first;
    int listenerCount = 0;
    
    /**
     * Class used to construct a two-way linked list of listeners.
     */
    public static class Linker implements java.io.Serializable {
        public Object listener;
        public Linker next, previous;
        public Linker(Object listener, Linker next, Linker previous) {
            this.listener = listener;
            this.next = next;
            this.previous = previous;
        }
    }

}
