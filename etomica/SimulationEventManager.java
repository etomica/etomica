package etomica;

/**
 * Class to take care of listener lists and event firing for simulation elements.
 * A class can make an instance of this manager as a field, and delegate any 
 * listener management functions to it.  This is done by the Controller class,
 * and could be used by Integrator.
 *
 * @see SimulationEvent
 * @see SimulationListener
 */
public class SimulationEventManager implements java.io.Serializable {
    
    //A Linker class is used to construct a linked list of listeners.
    private SimulationListener.Linker first, last;
    
    /**
     * Adds a listener.  Synchronized to avoid conflict with removeListener.
     */
    public synchronized void addListener(SimulationListener listener) {
        if(first == null) { //adding first listener
            first = new SimulationListener.Linker(listener, null, null);
            last = first;
        }
        else { //add listener to end of list
            last.next = new SimulationListener.Linker(listener, null, last);
            last = last.next;
        }
    }

    /**
     * Removes a listener.  Synchronized to avoid conflict with addListener.
     */
    public synchronized void removeListener(SimulationListener listener) {
        for(SimulationListener.Linker link=first; link!=null; link=link.next) {
            if(link.listener == listener) {
                if(link == first) {first = link.next;}
                if(link == last) {last = link.previous;}
                if(link.previous != null) link.previous.next = link.next;
                if(link.next != null) link.next.previous = link.previous;
                return;
            }
        }
    }

//probably ok to make not synchronized.  Addition and removal of listeners
//while firing event does not disturb complete event firing, since the the link keeps its next/previous
//handles when removed, and additions are made only to the end of the list
/*    public void fireEvent(SimulationEvent event) {
        for(SimulationEventListener.Linker link=first; link!=null; link=link.next) {
            link.listener.simulationAction(event);
        }
    }
*/
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

}
