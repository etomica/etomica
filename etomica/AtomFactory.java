package etomica;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 *
 * @author David Kofke
 */
public abstract class AtomFactory {
    
    protected final AtomReservoir reservoir;
    public final Simulation parentSimulation;
    private Species species;
    protected Configuration configuration;
    protected BondInitializer bondInitializer = BondInitializer.NULL;
    private Atom.AgentSource[] agentSource = new Atom.AgentSource[0];
    
    /**
     * @param sim the parent simulation using this factory
     * @param species the parent species using this factory (may be null)
     */
    public AtomFactory(Simulation sim) {
        parentSimulation = sim;
        reservoir = new AtomReservoir(this);
    }
    
    public Atom makeAtom() {
        return makeAtom((AtomTreeNodeGroup)reservoir.node);
    }
    
    public Atom makeAtom(AtomTreeNodeGroup parent) {
        Atom atom = reservoir.getAtom();
        if(atom == null) atom = build(parent);
        atom.node.setParent(parent);
        
        //add agents from any registered sources
        if(agentSource.length > 0) atom.agents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            atom.agents[i] = agentSource[i].makeAgent(atom);
        }
        
        return atom;
    }
    
    protected abstract Atom build(AtomTreeNodeGroup parent);
    
//    protected abstract void renew(Atom a);
    
//    public abstract boolean vetoAddition(Atom a); //be sure to check that a is non-null
    
    public AtomReservoir reservoir() {return reservoir;}
    
    public Simulation parentSimulation() {return parentSimulation;}
    
    public void setSpecies(Species species) {this.species = species;}
    public Species species() {return species;}
        
    public void setConfiguration(Configuration config) {configuration = config;}
    public Configuration getConfiguration() {return configuration;}
    
    public void setBondInitializer(BondInitializer bonder) {bondInitializer = bonder;}
    public BondInitializer getBondInitializer() {return bondInitializer;}
    
    /**
     * Adds given agent source to agent-source array and returns index
     * indicating where in atom agent-array the source's agent will
     * be placed.
     */
    public int requestAgentIndex(Atom.AgentSource aSource) {
        Atom.AgentSource[] newSource = new Atom.AgentSource[agentSource.length+1];
        for(int i=0; i<agentSource.length; i++) newSource[i] = agentSource[i];
        int index = agentSource.length;
        newSource[index] = aSource;
        agentSource = newSource;
        return index;
    }
}//end of AtomFactory