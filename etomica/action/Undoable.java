/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;
import etomica.Action;
/**
 * Interface for an action that can be "undone", such that its effect
 * is reversed if the undo() method is invoked immediately after the
 * actionPerformed().
 */
public interface Undoable extends Action {
	
	/**
	 * Causes the effect of most recent call to actionPerformed to be reversed.
	 * Assumes that no changes have been made in the interim that might
	 * alter the ability of the undo to recover the state before actionPerformed.
	 */
    public void undo();
}