package etomica; 

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
  
    public static String VERSION = "Potential1:01.07.26/"+Potential.VERSION;
    
    protected AtomIterator iterator;
    private Species species;
    
    public Potential1(PotentialGroup parent) {
        super(parent);
        iterator = new AtomIteratorSequential(true);
    }
    
    
    public Potential set(Atom a) {iterator.setBasis(a); return this;}

    public Potential set(Atom[] atoms) {
        if(atoms.length != 1) throw new IllegalArgumentException("Too many atoms in Potential1.set");
        return set(atoms[0]);
    }
    public Potential set(SpeciesMaster s) {
        if(species != null) iterator.setBasis(species.getAgent(s));
        else iterator.setBasis(s);
        return this;
      //  iterator.setBasis((species != null) ? species.getAgent(s.parentPhase()) : s);
    }
    
    public void setSpecies(Species s) {
        species = s;
    }
    
    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
                      
    public void setIterator(AtomIterator iterator) {
        this.iterator = iterator;
    }
    public AtomIterator iterator() {return iterator;}
        
}//end of Potential1



