package etomica.action.activity;

import etomica.action.Action;
import etomica.util.EnumeratedType;

public class ControllerEvent implements java.io.Serializable {
    
    protected final Controller controller;
    protected final Type type;
    protected final Action action;
    
    public ControllerEvent(Controller source, Type type) {
        this(source, type, null);
    }
    public ControllerEvent(Controller source, Type type, Action action) {
        this.controller = source;
        this.type = type;
        this.action = action;
    }
    
    public Action getAction() {return action;} 
    public Type getType() {return type;}
    public Controller getController() {return controller;}
    
    public static final Type START = new Type("Start",1);
    public static final Type START_ACTION = new Type("Start action",2);
    public static final Type END_ACTION = new Type("End action",3);
    public static final Type START_URGENT_ACTION = new Type("Start urgent action",4);
    public static final Type END_URGENT_ACTION = new Type("End urgent action",5);
    public static final Type NO_MORE_ACTIONS = new Type("No more actions",6);
    public static final Type HALTED = new Type("Halted",7);
    public static final Type RESET = new Type("Reset",8);
    
    public static class Type extends EnumeratedType {
        protected Type(String label, int index) {
            super(label);
            this.index = index;
        }       
        public static Type[] choices() { 
            return new Type[] {START,START_ACTION,END_ACTION,START_URGENT_ACTION,
                    END_URGENT_ACTION,NO_MORE_ACTIONS,HALTED,RESET};
        }
        public final int index;
    }//end of ValueType
}//end of ControllerEvent
