/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action.activity;

import etomica.action.IAction;
import etomica.util.EnumeratedType;
import etomica.util.IEvent;

public class ControllerEvent implements IEvent {
    
    protected final Controller controller;
    protected final Type type;
    protected final IAction action;
    
    public ControllerEvent(Controller source, Type type) {
        this(source, type, null);
    }
    public ControllerEvent(Controller source, Type type, IAction action) {
        this.controller = source;
        this.type = type;
        this.action = action;
    }
    
    public IAction getAction() {return action;} 
    public Type getType() {return type;}
    public Controller getController() {return controller;}

    public enum Type {
        START,
        START_ACTION,
        END_ACTION,
        START_URGENT_ACTION,
        END_URGENT_ACTION,
        NO_MORE_ACTIONS,
        HALTED,
        RESET
    }
}
