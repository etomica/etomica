package etomica;

/**
 * Superclass for all Potential classes.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 8/14/02 (DAK) made parentPotential mutable, so that potential can be added/removed from a potential group
  *               added setParentPotential for this purpose.
  */
public abstract class Potential extends SimulationElement {
    
    public static String VERSION = "Potential:01.07.22";
    
    private PotentialGroup parentPotential;
	private final PotentialAgent.List agentList = new PotentialAgent.List();
	protected final PotentialAgent.Iterator agentIterator = agentList.iterator();

    public Potential(PotentialGroup parent) {
        super(parent.parentSimulation(), Potential.class);
        parentPotential = parent;
        parentPotential.addPotential(this);
    }
    
    protected abstract PotentialAgent makeAgent();
    
    public void discardAgent(PotentialAgent agent) {
    	agentList.remove(agent);
    }   
    public PotentialAgent requestAgent() {
    	PotentialAgent agent = makeAgent();
    	agentList.add(agent);
    	return agent;
    }
    
    public final PotentialGroup parentPotential() {return parentPotential;}
    
    /**
     * Adds this potential to the given potential group.  No action is taken
     * if the new parent is the same as the current parent.  Otherwise, if
     * current parent is not null, this is removed from it.  Then this is
     * added via addPotential to the new parent.  If new parent is null, this
     * is removed from current parent and no new parent is set.
     */
    public void setParentPotential(PotentialGroup newParent) {
        if(newParent == parentPotential) return;
        if(parentPotential != null) parentPotential.removePotential(this);
        parentPotential = newParent;
        if(newParent != null) parentPotential.addPotential(this);
    }
    
    public abstract void setSpecies(Species[] s);
    public abstract Species[] getSpecies();
    /**
     * Marker interface for Null potentials, which are defined to have no action.
     */
    public interface Null {}
    
	final PotentialAgent.List agents = new PotentialAgent.List();
    
    
}//end of Potential