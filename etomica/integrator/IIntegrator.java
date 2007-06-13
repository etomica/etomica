package etomica.integrator;

import etomica.action.Action;
import etomica.exception.ConfigurationOverlapException;

public interface IIntegrator {

    /**
     * Performs the elementary integration step, such as a molecular dynamics
     * time step, or a Monte Carlo trial.
     */
    public void doStep();

    /**
     * Returns the number of steps performed by the integrator since it was
     * initialized.
     */
    public int getStepCount();

    /**
     * Defines the actions taken by the integrator to reset itself, such as
     * required if a perturbation is applied to the simulated phase (e.g.,
     * addition or deletion of a molecule). Also invoked when the
     * integrator is started or initialized.
     */
    //This should be called by subclasses /after/ they have performed their own
    //reset
    public void reset() throws ConfigurationOverlapException;

    /**
     * Initializes the integrator, performing the following steps: (1) deploys
     * agents in all atoms; (2) call reset method; (3) fires an event
     * indicating to registered listeners indicating that initialization has
     * been performed (i.e. fires IntervalEvent of type field set to
     * INITIALIZE).
     */
    public void initialize() throws ConfigurationOverlapException;

    /**
     * Returns true if initialize method has been called.
     */
    public boolean isInitialized();

    /**
     * Adds the given interval action to those that receive interval events
     * fired by this integrator.  Do not add an interval action that is already
     * an interval listener to this integrator.  Do not add a null listener.
     */
    public void addIntervalAction(Action newIntervalAction);

    /**
     * Sets the "priority" of the given intervalAction (which must have already
     * been added).  The priority determines the order in which interval
     * actions are performed.  Low priorities are fired first and high
     * priorities are fired last.  Priorities can range from 0 to 200.  If PBC
     * are enforced, they are enforced with a priority of 100.  Interval
     * actions with the same priority may be fired in any order.
     */
    public void setIntervalActionPriority(Action intervalAction, int newPriority);

    /**
     * Sets the number of integration steps between calls to the
     * intervalAction's actionPerformed method.  The interval must be positive.
     */
    public void setActionInterval(Action intervalAction, int newInterval);

    /**
     * Removes given interval action from those notified of interval events
     * fired by this integrator.  Do not try to remove an interval action that
     * is not a listener to this integrator.
     */
    public void removeIntervalAction(Action intervalAction);

    /**
     * Adds the given listener to those that receive non-interval events fired
     * by this integrator.  Do not add a listener that has been added or add a
     * null listener.
     */
    public void addNonintervalListener(IntegratorNonintervalListener iil);

    /**
     * Removes given listener from those notified of non-interval events fired
     * by this integrator.  Do not attempt to remove a listener which is not
     * listening to events from this integrator.
     */
    public void removeNonintervalListener(IntegratorNonintervalListener iil);

    /**
     * Returns all of the interval actions as an array.
     */
    public Action[] getIntervalActions();

    /**
     * Returns all of the non-interval listeners as an array.
     */
    public IntegratorNonintervalListener[] getNonintervalListeners();

}