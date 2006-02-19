package etomica.util;
import java.io.IOException;

import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationListener;

/**
 * Class to take care of listener lists and event firing for simulation elements.
 * A class can make an instance of this manager as a field, and delegate any 
 * listener management functions to it.  
 */
public abstract class EventManager implements java.io.Serializable {
    
    /**
     * Adds a listener.  Synchronized to avoid conflict with removeListener.
     */
    public synchronized void addListener(Object listener) {
        addListener(listener, true);
    }
    
    public synchronized void addListener(Object listener, boolean doSerialize) {
        if (listener.getClass().isInstance(getListenerClass())) {
            throw new IllegalArgumentException("must add listeners of class "+getListenerClass());
        }
        //add listener to beginning of list 
        //placement at end causes problem if a listener removes and then adds itself to the list as part of its response to the event
        first = new EventManager.Linker(listener, first, doSerialize);
    }

    /**
     * Removes a listener.  Synchronized to avoid conflict with addListener.
     */
    public synchronized void removeListener(Object listener) {
        Linker previous = null;
        for(Linker link=first; link!=null; link=link.next) {
            if(link.listener == listener) {
                if(link == first) {first = link.next;}
                if(previous != null) previous.next = link.next;
                return;
            }
            previous = link;
        }
    }
    
    protected abstract Class getListenerClass();

    protected transient EventManager.Linker first;
    
    private void writeObject(java.io.ObjectOutputStream out)
    throws IOException
    {
        // do nothing
        out.defaultWriteObject();
        
        // write # of listeners that will be serialized
        int count = 0;
        for (Linker link = first; link != null; link = link.next) {
            //look for non-transient listeners
            if (link.doSerialize)
                count++;
        }
        out.writeInt(count);

        for (Linker link = first; link != null; link = link.next) {
            //skip transient listeners
            if (link.doSerialize)
                out.writeObject(link.listener);
        }
    }

    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        // do nothing
        in.defaultReadObject();
        // we read the listener count but leave listenerCount at zero so we 
        // can invoke addListener
        int count = in.readInt();

        for (int i=0; i<count; i++) {
            addListener(in.readObject());
        }
    }
    
    /**
     * Class used to construct a two-way linked list of listeners.
     */
    protected static class Linker {
        public Object listener;
        public Linker next;
        protected boolean doSerialize;
        public Linker(Object listener, Linker next, boolean doSerialize) {
            this.listener = listener;
            this.next = next;
            this.doSerialize = doSerialize;
        }
    }

}
