package etomica.integrator;

import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.api.IAction;
import etomica.api.IIntegrator;
import etomica.api.IVectorMutable;

/**
 * Integrator implements the algorithm used to move the atoms around and
 * generate new configurations in one or more boxs. All integrator techniques,
 * such as molecular dynamics or Monte Carlo, are implemented via subclasses of
 * this Integrator class. The Integrator's activities are managed via the
 * actions of the governing Controller.
 * 
 * @author David Kofke and Andrew Schultz
 */
public abstract class Integrator implements java.io.Serializable, IIntegrator {

    private static final long serialVersionUID = 1L;
    protected boolean initialized = false;
    private final LinkedList<ListenerWrapper> intervalListeners = new LinkedList<ListenerWrapper>();
    private ListenerWrapper[] listenerWrapperArray = new ListenerWrapper[0];
    protected int interval;
    private int iieCount;
    protected long stepCount;

    public Integrator() {
        setEventInterval(1);
        stepCount = 0;
    }

    /**
     * @return value of interval (number of doStep calls) between
     * firing of interval events by integrator.
     */
    public int getEventInterval() {
        return interval;
    }
    
    /**
     * Sets value of interval between successive firing of integrator interval events.
     * @param interval
     */
    public void setEventInterval(int interval) {
        this.interval = interval;
        if (iieCount > interval) {
            iieCount = interval;
        }
    }
    
    public final void doStep() {
        stepCount++;
        doStepInternal();
        if(--iieCount == 0) {
            fireIntervalEvent();
            iieCount = interval;
        }
    }
    
    /**
     * Performs the elementary integration step, such as a molecular dynamics
     * time step, or a Monte Carlo trial.
     */
    protected abstract void doStepInternal();

    /**
     * Returns the number of steps performed by the integrator since it was
     * initialized.
     */
    public long getStepCount() {
        return stepCount;
    }
    
    /**
     * Defines the actions taken by the integrator to reset itself, such as
     * required if a perturbation is applied to the simulated box (e.g.,
     * addition or deletion of a molecule). Also invoked when the
     * integrator is started or initialized.
     */
    //This should be called by subclasses before they have performed their own
    //reset
    public void reset() {
        if (!initialized) {
            setup();
        }
    }

    public void resetStepCount() {
        stepCount = 0;
    }

    /**
     * Returns true if initialize method has been called.
     */
    public boolean isInitialized() {
        return initialized;
    }

    /**
     * Perform initialization.  Subclasses can override this method to set up
     * before integration begins.
     */
    protected void setup() {
        resetStepCount();
        iieCount = interval;
        initialized = true;
    }

    /**
     * Arranges interval listeners registered with this iterator in order such that
     * those with the smallest (closer to zero) priority value are performed
     * before those with a larger priority value.  This is invoked automatically
     * whenever a listener is added or removed.
     */
    protected synchronized void sortListeners() {
        //sort using linked list, but put into array afterwards
        //for rapid looping (avoid repeated construction of iterator)
        Collections.sort(intervalListeners);
        listenerWrapperArray = intervalListeners.toArray(new ListenerWrapper[0]);
    }
    
    /**
     * Removes all interval and noninterval listeners.
     */
    public synchronized void removeAllListeners() {
        intervalListeners.clear();
        sortListeners();
    }
    
    /**
     * Notifies registered listeners that an interval has passed. Not
     * synchronized, so unpredictable behavior if listeners are added while
     * notification is in process (this should be rare).
     */
    protected void fireIntervalEvent() {
        for(int i=0; i<listenerWrapperArray.length; i++) {
            if (--listenerWrapperArray[i].intervalCount == 0) {
                listenerWrapperArray[i].intervalAction.actionPerformed();
                listenerWrapperArray[i].intervalCount = listenerWrapperArray[i].interval;
            }
        }
    }

    /**
     * Adds the given interval listener to those that receive interval events fired by
     * fired by this integrator.  If the action has already been added to
     * integrator, an exception is thrown.  If listener is null,
     * NullPointerException is thrown.
     */
    public synchronized void addIntervalAction(IAction newIntervalAction) {
        if(newIntervalAction == null) throw new NullPointerException("Cannot add null as a listener to Integrator");
        if (intervalListeners.contains(newIntervalAction)) {
            throw new RuntimeException(newIntervalAction+" is already an interval action");
        }
        intervalListeners.add(new ListenerWrapper(newIntervalAction, 100));
        sortListeners();
    }
        
    /**
     * Sets the priority for the given interval action.  The order in which
     * interval actions for a given time step are fired is based on the
     * priority.  Actions with a small priority (closer to zero) are invoked
     * before those with a larger priority.  The priority must be positive and
     * not greater than 200.  The default priority is 100.
     */
    public synchronized void setIntervalActionPriority(IAction intervalAction, int newPriority) {
        if (newPriority < 0 || newPriority > 200) {
            throw new RuntimeException("Priority may not be negative or greater than 200");
        }
        findWrapper(intervalAction).priority = newPriority;
        sortListeners();
    }
        
    /**
     * Sets the interval for the given action.  After this call, the action
     * will be fired after every N steps, where N is the newInterval given to
     * this method.
     */
    public synchronized void setActionInterval(IAction intervalAction, int newInterval) {
        if (newInterval < 1) {
            throw new RuntimeException("Action interval must be positive");
        }
        findWrapper(intervalAction).setInterval(newInterval);
    }
        
    /**
     * Returns the interval for the given action.  The action is fired after
     * every N steps, where N is the number returned by this method.
     */
    public synchronized int getActionInterval(IAction intervalAction) {
        return findWrapper(intervalAction).interval;
    }
        
    /**
     * Finds and returns the ListenerWrapper used to put the given listener in the list.
     * Returns null if listener is not in list.
     */
    private ListenerWrapper findWrapper(IAction intervalAction) {
        Iterator<ListenerWrapper> iterator = intervalListeners.iterator();
        while(iterator.hasNext()) {
            ListenerWrapper wrapper = iterator.next();
            if(wrapper.intervalAction == intervalAction) return wrapper;//found it
        }
        return null;//didn't find it in list      
    }

    /**
     * Removes given interval listener from those notified of interval events fired
     * by this integrator.  No action results if given listener is null or is not registered
     * with this integrator.
     */
    public synchronized void removeIntervalAction(IAction intervalAction) {
        intervalListeners.remove(findWrapper(intervalAction));
        sortListeners();
    }

    public synchronized int getIntervalActionCount() {
    	return listenerWrapperArray.length;
    }

    public synchronized IAction getIntervalAction(int index) {
        return listenerWrapperArray[index].intervalAction;
    }

    /**
     * Wraps interval listener an implements a Comparable interface based
     * on the listeners priority value.
     * This class has a natural ordering that is inconsistent with equals.
     */
    private static class ListenerWrapper implements Comparable<ListenerWrapper>, java.io.Serializable {
        private static final long serialVersionUID = 1L;
        protected final IAction intervalAction;
        protected int priority;
        protected int interval;
        protected int intervalCount;
        
        protected ListenerWrapper(IAction intervalAction, int priority) {
            this.intervalAction = intervalAction;
            this.priority = priority;
            interval = 1;
            intervalCount = 1;
        }
        
        public void setInterval(int newInterval) {
            interval = newInterval;
            intervalCount = interval;
        }
        
        public int compareTo(ListenerWrapper obj) {
            int objPriority = obj.priority;
            //we do explicit comparison of values (rather than returning
            //their difference) to avoid potential problems with large integers.
            if(priority < objPriority) return -1;
            if(priority == objPriority) return 0;
            return +1;
         }
    }

    /**
     * Integrator agent that holds a force vector. Used to indicate that an atom
     * could be under the influence of a force.
     */
    public interface Forcible {
        public IVectorMutable force();
    }

    /**
     * Integrator agent that holds a torque vector.
     */
    public interface Torquable {
        public IVectorMutable torque();
    }
}