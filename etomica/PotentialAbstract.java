package etomica;

/**
 * Superclass for all Potential classes
 */
public abstract class PotentialAbstract implements Simulation.Element, java.io.Serializable {
    
    public static String VERSION = "PotentialAbstract:01.06.12";
    
    private final Simulation parentSimulation;
    private boolean added = false;
    private String name, label;
    protected PotentialAbstract next;
//    private PotentialAbstract previous;
    
    public PotentialAbstract(Simulation sim) {
        parentSimulation = sim;
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return PotentialAbstract.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    public final String getName() {return name;}
    public final void setName(String name) {this.name = name;}
    
    public final String getLabel() {return label;}
    public final void setLabel(String text) {label = text;}
    public String toString() {return label;}
    
    public abstract double energy();
    public abstract double energy(Atom atom);
//    public abstract double energy(AtomGroup group);
    
    public PotentialAbstract next() {return next;}
//    public PotentialAbstract previous() {return previous;}
    public void setNext(PotentialAbstract potl) {next = potl;}
    
    public interface Hard {
        public void findCollisions(IntegratorHardAbstract.CollisionHandler c);
        public void findCollisions(Atom a, IteratorDirective id, IntegratorHardAbstract.CollisionHandler c);
    }
}