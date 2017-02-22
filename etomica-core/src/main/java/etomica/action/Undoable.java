/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;


/**
 * Interface for an action that can be "undone", such that its effect
 * is reversed if the undo() method is invoked immediately after the
 * actionPerformed().
 */
public interface Undoable extends IAction {
	
	/**
	 * Causes the effect of most recent call to actionPerformed to be reversed.
	 * Assumes that no changes have been made in the interim that might
	 * alter the ability of the undo to recover the state before actionPerformed.
	 */
    public void undo();
}