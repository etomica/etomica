package etomica;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 *
 * @author David Kofke
 */
public abstract class AtomFactory {
    
    protected final AtomReservoir reservoir;
    protected Simulation parentSimulation;
    protected Configuration configuration;
    protected BondInitializer bondInitializer = BondInitializer.NULL;
    
    public AtomFactory(Simulation sim) {
        parentSimulation = sim;
        reservoir = new AtomReservoir(sim);
    }
    
    public Atom makeAtom() {
        Atom atom = reservoir.removeAtom();
        if(atom == null) atom = build();
        return atom;
    }
    
    protected abstract Atom build();
    
    protected abstract void renew(Atom a);
    
    public abstract boolean vetoAddition(Atom a); //be sure to check that a is non-null
    
    public AtomReservoir reservoir() {return reservoir;}
    
    public void setConfiguration(Configuration config) {configuration = config;}
    public Configuration getConfiguration() {return configuration;}
    
    public void setBondInitializer(BondInitializer bonder) {bondInitializer = bonder;}
    public BondInitializer getBondInitializer() {return bondInitializer;}
    
}