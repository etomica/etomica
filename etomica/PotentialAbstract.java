package etomica;

//Java2 imports
//import java.util.HashMap;

import etomica.utility.HashMap;

/**
 * Superclass for all Potential classes
 */
public abstract class PotentialAbstract implements Simulation.Element, java.io.Serializable {
    
    public static String VERSION = "PotentialAbstract:01.06.28";
    
    private final Simulation parentSimulation;
    private boolean added = false;
    private String name, label;
    
    /**
     * Hashtable to associate agents with phases
     * Agent is placed in hashmap of potential in the constructor of PotentialAgent.
     */
    final HashMap agents = new HashMap();

    public PotentialAbstract(Simulation sim) {
        parentSimulation = sim;
        parentSimulation.register(this);
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

    /**
     * Returns the agent of this species in the given phase.
     * Hashmap is used to connect phase(key)-agent(value) pairs.
     * Agent is placed in hashmap of potential in the agent's constructor.
     * 
     * @param p The phase for which this species' agent is requested
     * @return The agent of this species in the phase
     */
    public final PotentialAgent getAgent(Phase p) {return (PotentialAgent)agents.get(p);}
    
    public final HashMap agents() {return agents;}
 
}//end of PotentialAbstract