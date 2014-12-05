/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.IAction;

/**
 * Wrapper of an etomica Action that permits it to be used as a java.awt action listener.
 * Useful, for example, to connect a button to the action.
 *
 * @author David Kofke
 */
public class ActionGraphic implements IAction, java.awt.event.ActionListener, java.io.Serializable { 
    
    public ActionGraphic(IAction action) {
        simulationAction = action;
    }
    
    public void actionPerformed(java.awt.event.ActionEvent evt) {
        simulationAction.actionPerformed();
    }
    
    public void actionPerformed() {
        simulationAction.actionPerformed();
    }
    
    private static final long serialVersionUID = 1L;
    private final IAction simulationAction;
}
