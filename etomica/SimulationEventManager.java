package etomica;
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
    
    //A Linker class is used to construct a linked list of listeners.
    private SimulationListener.Linker first;
    int listenerCount = 0;
    
    /**
     * Adds a listener.  Synchronized to avoid conflict with removeListener.
     */
    public synchronized void addListener(SimulationListener listener) {
        //add listener to beginning of list 
        //placement at end causes problem if a listener removes and then adds itself to the list as part of its response to the event
        first = new SimulationListener.Linker(listener, first, null);
        if(first.next != null) first.next.previous = first;
        listenerCount++;
    }

    /**
     * Removes a listener.  Synchronized to avoid conflict with addListener.
     */
    public synchronized void removeListener(SimulationListener listener) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
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
     public SimulationListener.Linker first() {return first;}

//probably ok to make not synchronized.  Addition and removal of listeners
//while firing event does not disturb complete event firing, since the the link keeps its next/previous
//handles when removed, and additions are made only to the end of the list
    public void fireEvent(SimulationEvent event) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
            link.listener.actionPerformed(event);
        }
    }

    public void fireEvent(PhaseEvent event) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
            ((PhaseListener)link.listener).actionPerformed(event);
        }
    }

    public void fireEvent(ControllerEvent event) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
            ((ControllerListener)link.listener).actionPerformed(event);
        }
    }

    public void fireEvent(MCMoveEvent event) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
            ((MCMoveListener)link.listener).actionPerformed(event);
        }
    }

    public void fireEvent(LatticeEvent event) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
            ((LatticeListener)link.listener).actionPerformed(event);
        }
    }

}
