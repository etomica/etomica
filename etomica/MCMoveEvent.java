package etomica;

public class MCMoveEvent extends SimulationEvent {
    
//    public final Type type;
    public MCMove mcMove;
    public boolean isTrialNotify;
    public boolean wasAccepted;
    
    public MCMoveEvent(Object source) {
        super(source);
    }
    
    /**
     * Enumerated type that indicate the type of MCMoveEvent.
     */
/*	public static class Type extends Constants.TypedConstant {
        private Type(String label) {super(label);}       
        public Constants.TypedConstant[] choices() {return (Constants.TypedConstant[])CHOICES;}
        public static final Type[] CHOICES = 
            new Type[] {
                new Type("Trial"),
                new Type("Accept"),
                new Type("Reject")};
    }//end of Type
    public static final Type TRIAL = Type.CHOICES[0];
    public static final Type ACCEPT = Type.CHOICES[1];
    public static final Type REJECT = Type.CHOICES[2];
    */
}//end of MCMoveEvent