package etomica.event;

import javax.swing.event.*;

/**
 * Class to manage listeners of change events and fire events to them.
 */
 
 /* History
  * 10/12/02 (DAK) new
  */
  
public class ChangeEventManager {
    
    EventListenerList listenerList = new EventListenerList();
    final ChangeEvent changeEvent;
    
    /**
     * Argument to constructor will be the source of events fired by the
     * manager if fireEvent is invoke with no argument.
     */
    public ChangeEventManager(Object source) {
        changeEvent = new ChangeEvent(source);
    }

    public void addChangeListener(ChangeListener l) {
        listenerList.add(ChangeListener.class, l);
    }

    public void removeChangeListener(ChangeListener l) {
        listenerList.remove(ChangeListener.class, l);
    }

     
    // Notify all listeners that have registered interest for
    // notification on this event type.  The event instance 
    // is lazily created using the parameters passed into 
    // the fire method.

    public void fireChangeEvent() {
        fireChangeEvent(changeEvent);
    }
    public void fireChangeEvent(ChangeEvent evt) {
        // Guaranteed to return a non-null array
        Object[] listeners = listenerList.getListenerList();
        // Process the listeners last to first, notifying
        // those that are interested in this event
        for (int i = listeners.length-2; i>=0; i-=2) {
	        if (listeners[i]==ChangeListener.class) {
	        ((ChangeListener)listeners[i+1]).stateChanged(evt);
	        }	       
        }
    }	
}//end of ChangeEventManager