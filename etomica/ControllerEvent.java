package etomica;

public class ControllerEvent extends SimulationEvent {
    
    protected Controller controller;
    protected Type type;
    
    public ControllerEvent(Object source) {
        this(source, null);
    }
    public ControllerEvent(Object source, Type t) {
        super(source);
        type = t;
    }
    
    public void setType(Type t) {type = t;}
    public Type type() {return type;}
    
    public final ControllerEvent setController(Controller c) {controller = c; return this;}
    public final Controller controller() {return controller;}
    
    //no Type classes yet defined
    public static class Type extends Constants.TypedConstant {
        private Type(String label) {super(label);}
        public static final Type[] CHOICES = new Type[] {};
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Type   
}//end of ControllerEvent