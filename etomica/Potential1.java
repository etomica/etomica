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
    
    public void setSpecies(Species[] species) {
        switch (species.length) {
            case 1: setSpecies(species[0]);
                    break;
            default: throw new IllegalArgumentException("Wrong number of species given in Potential2");
        }
    }
    /**
     * Returns an array of length 2 with the species to which this potential applies.
     * Returns null if no species has been set, which is the case if the potential
     * is not describing interactions between molecule-level Atoms.
     */
    public Species[] getSpecies() {
        if(species == null) return null;
        else return new Species[] {species};
    }

    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
                      
    public void setIterator(AtomIterator iterator) {
        this.iterator = iterator;
    }
    public AtomIterator iterator() {return iterator;}
    
    /**
     * Marker interface indicating that a one-body potential is an intramolecular
     * potential, and not, e.g., a potential of interaction with an external field.
     * This is useful when computing energy changes for molecule translations and
     * rotations, for which intramolecular contributions can be ignored.
     */
    public interface Intramolecular {}
        
}//end of Potential1



