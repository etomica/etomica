package etomica;

public abstract class BondInitializer {
    
    public abstract void makeBonds(Atom a);
    
    
    /**
     * Initializer that causes no bonds to be made.
     */
    public static final BondInitializer NULL = new Null();
    private static final class Null extends BondInitializer {
        public void makeBonds(Atom a) {}
    }
}//end of BondInitializer