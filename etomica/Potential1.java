package etomica; 

/**
 * Potential acting on a single atom.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
  
    public static String VERSION = "Potential1:01.06.27/"+Potential.VERSION;
    
    protected AtomIterator iterator;
    
    public Potential1(Simulation sim) {
        super(sim);
    }
    
    public Potential set(Atom a) {iterator.setBasis(a); return this;}

    public Potential set(Atom a1, Atom a2) {return null;} //throw an exception? redesign?
    
    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
          
    /**
     * Default iterator gives no atoms.
     */
    protected void makeDefaultIterator() {
        iterator = AtomIterator.NULL;
    }
            
    public void setIterator(AtomIterator iterator) {
        this.iterator = iterator;
    }
    public AtomIterator iterator() {return iterator;}
        
}//end of Potential1



