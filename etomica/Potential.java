package etomica;

/**
 * Superclass for all Potential classes
 *
 * @author David Kofke
 */
public abstract class Potential implements Simulation.Element, java.io.Serializable {
    
    public static String VERSION = "Potential:01.07.22";
    
    private final Simulation parentSimulation;
    private boolean added = false;
    private String name, label;
    
    public Potential(Simulation sim) {
        parentSimulation = sim;
        parentSimulation.register(this);
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Potential.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    public final String getName() {return name;}
    public final void setName(String name) {this.name = name;}
    
    public final String getLabel() {return label;}
    public final void setLabel(String text) {label = text;}
    public String toString() {return label;}
        
    public abstract void calculate(IteratorDirective id, PotentialCalculation pc);
    
    //Sets the basis for iteration
    public abstract Potential set(Atom a);
    public abstract Potential set(Atom a1, Atom a2);
    
}//end of Potential