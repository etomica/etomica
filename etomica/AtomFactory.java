package etomica;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 *
 * @author David Kofke
 */
public abstract class AtomFactory {
    
    private final AtomReservoir reservoir;
    protected Simulation parentSimulation;
    
    public AtomFactory(Simulation sim) {
        parentSimulation = sim;
        reservoir = new AtomReservoir(sim);
    }
    
    public Atom makeAtom() {
        Atom atom = reservoir.removeAtom();
        if(atom == null) atom = build(reservoir);
        return atom;
    }
    
    protected abstract Atom build(AtomGroup parent);
    
    protected abstract void renew(Atom a);
    
    public abstract boolean vetoAddition(Atom a); //be sure to check that a is non-null
    
    public AtomReservoir reservoir() {return reservoir;}
    
}