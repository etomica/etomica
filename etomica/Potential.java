package etomica;

/**
 * Superclass for all Potential classes
 *
 * @author David Kofke
 */
public abstract class Potential extends SimulationElement {
    
    public static String VERSION = "Potential:01.07.22";
    
    private final PotentialGroup parentPotential;
    
    public Potential(PotentialGroup parent) {
        super(parent.parentSimulation(), Potential.class);
        parentPotential = parent;
        parentPotential.addPotential(this);
    }
        
    public abstract void calculate(IteratorDirective id, PotentialCalculation pc);
    
    //Sets the basis for iteration
    public abstract Potential set(Atom a);
    public abstract Potential set(Atom a1, Atom a2);
    public abstract Potential set(SpeciesMaster s);    
    /**
     * Marker interface for Null potentials, which are defined to have no action.
     */
    public interface Null {}
    
}//end of Potential