package etomica; 

/**
 * Potential acting on a pair of atoms or atom groups.
 *
 * @author David Kofke
 */
public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:01.07.03/"+Potential.VERSION;
    
    protected AtomPairIterator iterator;
    protected Space.Vector work1;
    
    public Potential2(Simulation sim) {
        super(sim);
        work1 = sim.space().makeVector();
    }
    
    public abstract double energy(AtomPair pair);

    public Potential set(Atom a) {return set(a,a);}
    public Potential set(Atom a1, Atom a2) {iterator.setBasis(a1, a2); return this;}
 
    /**
     * Default iterator yields no pairs.
     */
    protected void makeDefaultIterator() {
        iterator = AtomPairIterator.NULL;
    }
            
    public void setIterator(AtomPairIterator iterator) {
        this.iterator = iterator;
    }
    public AtomPairIterator iterator() {return iterator;}
    
}//end of Potential2



