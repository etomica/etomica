package etomica;

/**
 * Class to take care of listener lists and event firing for simulation elements.
 * A class can make an instance of this manager as a field, and delegate any 
 * listener management functions to it.  This is done by the Controller class,
 * and could be used by Integrator.
 *
 * @see SimulationEvent
 * @see SimulationEventListener
 */
public class SimulationEventManager implements java.io.Serializable {
    
    //A Linker class is used to construct a linked list of listeners. (say that 10 times quickly!)
    private SimulationEventListener.Linker first, last;
    
    /**
     * Adds a listener.  Synchronized to avoid conflict with removeListener.
     */
    public synchronized void addListener(SimulationEventListener listener) {
        if(first == null) { //adding first listener
            first = new SimulationEventListener.Linker(listener, null, null);
            last = first;
        }
        else { //add listener to end of list
            last.next = new SimulationEventListener.Linker(listener, null, last);
            last = last.next;
        }
    }

    /**
     * Removes a listener.  Synchronized to avoid conflict with addListener.
     */
    public synchronized void removeListener(SimulationEventListener listener) {
        for(SimulationEventListener.Linker link=first; link!=null; link=link.next) {
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
    public void fireEvent(SimulationEvent event) {
        for(SimulationEventListener.Linker link=first; link!=null; link=link.next) {
            link.listener.simulationAction(event);
        }
    }
}
