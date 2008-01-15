package etomica.integrator;

import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.action.Action;
import etomica.exception.ConfigurationOverlapException;
import etomica.space.IVector;

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

    protected boolean initialized = false;
    private final LinkedList intervalListeners = new LinkedList();
    private final LinkedList listeners = new LinkedList();
    private ListenerWrapper[] listenerWrapperArray = new ListenerWrapper[0];
    protected int interval;
    private int iieCount;
    protected int stepCount;

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
    public int getStepCount() {
        return stepCount;
    }
    
    /**
     * Defines the actions taken by the integrator to reset itself, such as
     * required if a perturbation is applied to the simulated box (e.g.,
     * addition or deletion of a molecule). Also invoked when the
     * integrator is started or initialized.
     */
    //This should be called by subclasses /after/ they have performed their own
    //reset
    public void reset() throws ConfigurationOverlapException {
        fireNonintervalEvent(new IntegratorNonintervalEvent(this, IntegratorNonintervalEvent.RESET));
    }

    /**
     * Initializes the integrator, performing the following steps: (1) deploys
     * agents in all atoms; (2) call reset method; (3) fires an event
     * indicating to registered listeners indicating that initialization has
     * been performed (i.e. fires IntervalEvent of type field set to
     * INITIALIZE).
     */
    public final void initialize() throws ConfigurationOverlapException {
        stepCount = 0;
        initialized = false;
        setup();
        initialized = true;
        reset();
        iieCount = interval;
    }
    
    /**
     * Returns true if initialize method has been called.
     */
    public boolean isInitialized() {
        return initialized;
    }

    protected void setup() throws ConfigurationOverlapException {}

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
        listenerWrapperArray = (ListenerWrapper[])intervalListeners.toArray(new ListenerWrapper[0]);
    }
    
    /**
     * Removes all interval and noninterval listeners.
     */
    public synchronized void removeAllListeners() {
        intervalListeners.clear();
        listeners.clear();
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
     * Notifies registered listeners of a non-interval event. 
     */
    protected synchronized void fireNonintervalEvent(IntegratorNonintervalEvent ie) {
        Iterator iter = listeners.iterator();
        while(iter.hasNext()) {
            ((IntegratorNonintervalListener)iter.next()).nonintervalAction(ie);
        }
    }

    /**
     * Adds the given interval listener to those that receive interval events fired by
     * fired by this integrator.  If the action has already been added to
     * integrator, an exception is thrown.  If listener is null,
     * NullPointerException is thrown.
     */
    public synchronized void addIntervalAction(Action newIntervalAction) {
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
    public synchronized void setIntervalActionPriority(Action intervalAction, int newPriority) {
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
    public synchronized void setActionInterval(Action intervalAction, int newInterval) {
        if (newInterval < 1) {
            throw new RuntimeException("Action interval must be positive");
        }
        findWrapper(intervalAction).setInterval(newInterval);
    }
        
    /**
     * Returns the interval for the given action.  The action is fired after
     * every N steps, where N is the number returned by this method.
     */
    public synchronized int getActionInterval(Action intervalAction) {
        return findWrapper(intervalAction).interval;
    }
        
    /**
     * Finds and returns the ListenerWrapper used to put the given listener in the list.
     * Returns null if listener is not in list.
     */
    private ListenerWrapper findWrapper(Action intervalAction) {
        Iterator iterator = intervalListeners.iterator();
        while(iterator.hasNext()) {
            ListenerWrapper wrapper = (ListenerWrapper)iterator.next();
            if(wrapper.intervalAction == intervalAction) return wrapper;//found it
        }
        return null;//didn't find it in list      
    }

    /**
     * Removes given interval listener from those notified of interval events fired
     * by this integrator.  No action results if given listener is null or is not registered
     * with this integrator.
     */
    public synchronized void removeIntervalAction(Action intervalAction) {
        intervalListeners.remove(findWrapper(intervalAction));
        sortListeners();
    }

    /**
     * Adds the given listener to those that receive non-interval events fired by
     * this integrator.  If listener has already been added to integrator, it is
     * not added again.  If listener is null, NullPointerException is thrown.
     */
    public synchronized void addNonintervalListener(IntegratorNonintervalListener iil) {
        if(iil == null) throw new NullPointerException("Cannot add null as a listener to Integrator");
        if(!listeners.contains(iil)) {
            listeners.add(iil);
        }
    }
    
    /**
     * Removes given listener from those notified of non-interval events fired
     * by this integrator.  No action results if given listener is null or is not registered
     * with this integrator.
     */
    public synchronized void removeNonintervalListener(IntegratorNonintervalListener iil) {
        listeners.remove(iil);
    }

    public synchronized Action[] getIntervalActions() {
        Action[] intervalListenerArray = new Action[listenerWrapperArray.length];
        for (int i=0; i<intervalListenerArray.length; i++) {
            intervalListenerArray[i] = listenerWrapperArray[i].intervalAction;
        }
        return intervalListenerArray;
    }

    public synchronized IntegratorNonintervalListener[] getNonintervalListeners() {
        IntegratorNonintervalListener[] listenerArray = new IntegratorNonintervalListener[listeners.size()];
        Iterator iter = listeners.iterator();
        int i=0;
        while(iter.hasNext()) {
            listenerArray[i++] = (IntegratorNonintervalListener)iter.next();
        }
        return listenerArray;
    }

    /**
     * Wraps interval listener an implements a Comparable interface based
     * on the listeners priority value.
     * This class has a natural ordering that is inconsistent with equals.
     */
    private static class ListenerWrapper implements Comparable, java.io.Serializable {
        private static final long serialVersionUID = 1L;
        protected final Action intervalAction;
        protected int priority;
        protected int interval;
        protected int intervalCount;
        
        protected ListenerWrapper(Action intervalAction, int priority) {
            this.intervalAction = intervalAction;
            this.priority = priority;
            interval = 1;
            intervalCount = 1;
        }
        
        public void setInterval(int newInterval) {
            interval = newInterval;
            intervalCount = interval;
        }
        
        public int compareTo(Object obj) {
            int objPriority = ((ListenerWrapper)obj).priority;
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
        public IVector force();
    }
}