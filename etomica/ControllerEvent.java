package etomica;

import etomica.Constants.TypedConstant;
import etomica.action.Action;

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
    
    public static class Type extends TypedConstant {
        protected Type(String label, int index) {
            super(label);
            this.index = index;
        }       
        public Constants.TypedConstant[] choices() {return CHOICES;}
        public final int index;
    }//end of ValueType
    protected static final Type[] CHOICES = 
        new Type[] {
            new Type("Start", 0),
            new Type("Start action", 1), 
            new Type("End action", 2),
            new Type("Start urgent action", 3),
            new Type("End urgent action", 4),
            new Type("No more actions", 5),
            new Type("Halted", 6)};
    public static final Type START = CHOICES[0];
    public static final Type START_ACTION = CHOICES[1];
    public static final Type END_ACTION = CHOICES[2];
    public static final Type START_URGENT_ACTION = CHOICES[3];
    public static final Type END_URGENT_ACTION = CHOICES[4];
    public static final Type NO_MORE_ACTIONS = CHOICES[5];
    public static final Type HALTED = CHOICES[6];
}//end of ControllerEvent
