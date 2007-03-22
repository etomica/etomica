package etomica.action;

/**
 * Interface for a class that performs a specific action. If the action is
 * performed quickly, the class may simply implement this interface. If the
 * action is expected to require a significant amount of time to complete, the
 * class should extend Activity, which implements this interface and provides
 * methods to support pausing/resuming/terminating the action.
 */
public interface Action {

    /**
     * Completes the action defined by the class implementing this interface.
     */
	public void actionPerformed();

}