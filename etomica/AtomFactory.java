package etomica;

import etomica.atom.AtomSequencerFactory;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 *
 * @author David Kofke
 */
 
public abstract class AtomFactory {
    
    public final Space space;
    protected Conformation conformation;
    protected AtomSequencerFactory sequencerFactory;
    protected AtomTreeNodeFactory nodeFactory;
    private Atom.AgentSource[] agentSource = new Atom.AgentSource[0];
    protected final AtomType atomType;
    
    /**
     * Makes an atom factory with atoms having AtomSequencerSimple and
     * AtomTreeNodeGroup for sequencer and node, respectively.
     * @param space
     */
    public AtomFactory(Space space, AtomType atomType) {
        this(space, atomType, AtomSequencerFactory.SIMPLE);
    }
    
    public AtomFactory(Space space, AtomType atomType, AtomSequencerFactory sequencerFactory) {
    	this(space, atomType, sequencerFactory, AtomTreeNodeGroup.FACTORY);
    }
    
    public AtomFactory(Space space, AtomType atomType, AtomSequencerFactory sequencerFactory, AtomTreeNodeFactory nodeFactory) {
        this.space = space;
        this.sequencerFactory = sequencerFactory;
        this.nodeFactory = nodeFactory;
        this.atomType = atomType;
        atomType.creator = this;
    }
    
    /**
     * Builds and returns the atom/atomgroup made by this factory.
     * Implementation of this method in the subclass defines the 
     * product of this factory.
     */
    public abstract Atom makeAtom();

    /**
     * Identifies the species for which this factory makes its atoms.
     * Should be invoked only in the species constructor, and by any
     * an atom factory on its child factories.
     * @param species
     */
    public abstract void setSpecies(Species species);
    
    /**
     * Returns the species that is using this factory or its parent factory.
     * @return
     */
    public Species getSpecies() {return atomType.getSpecies();}
        
    /**
     * Method used by subclasses to make the root atom of the group it is building.
     */
    protected Atom newParentAtom() {
        Atom atom = new Atom(space, atomType, nodeFactory, sequencerFactory);
        //add agents from any registered sources
        if(agentSource.length > 0) atom.agents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            atom.agents[i] = agentSource[i].makeAgent(atom);
        }
        return atom;
    }

    /**
     * Returns the atomType instance given to all atoms made by this factory.
     */
    public AtomType getType() {
        return atomType;
    }
    
    /**
     * Returns the space used to build the atoms made by this factory.
     */
    public Space getSpace() {
        return space;
    }

    /**
     * Sets the conformation used to set the standard arrangement of
     * the atoms/atom-groups produced by this factory.
     */
    public void setConformation(Conformation config) {
        conformation = config;
    }
    
    /**
     * Returns the conformation used to set the standard arrangement of
     * the atoms/atom-groups produced by this factory.
     */
    public Conformation getConformation() {return conformation;}
    
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