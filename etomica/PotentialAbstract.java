package etomica;

/**
 * Superclass for all Potential classes
 */
public abstract class PotentialAbstract implements Simulation.Element, java.io.Serializable {
    
    public static String VERSION = "PotentialAbstract:01.06.24";
    
    private final Simulation parentSimulation;
    private boolean added = false;
    private String name, label;
    
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
    
    public abstract PotentialAgent makeAgent(Phase p);
    
/*    public abstract void deployAgent(Phase phase) {
        PotentialAgent agent = makeAgent(phase);
        if(phase.potential == null) {
            phase.potential = agent;
        }
        else if(phase.potential instanceof PotentialGroup) {
            phase.potential
 */   
}//end of PotentialAbstract