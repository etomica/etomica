package etomica; 

/**
 * Potential acting on a trio of atoms or atom groups.
 *
 * @author David Kofke
 */
public abstract class Potential3 extends Potential {
  
    public static String VERSION = "Potential2:01.08.06/"+Potential.VERSION;
    
    protected Atom3Iterator iterator;
    private Species species1, species2, species3;
    
    public Potential3(PotentialGroup parent) {
        super(parent);
//        iterator = new Atom3Iterator(parentSimulation().space());
    }
    
    public abstract double energy(Atom3 atom3);
    
    public void setSpecies(Species s1, Species s2, Species s3) {//throw exception if either are null
		setSpecies(new Species[] {s1, s2, s3});
    }

    public void setSpecies(Species[] species) {
		species1 = species[0];
		species2 = species[1];
		species3 = species[2];
        switch (species.length) {
            case 1: 
                    break;
            case 2:
            		break;
            case 3: 
                    break;
            default: throw new IllegalArgumentException("Wrong number of species given in Potential3");
        }
		if(species1 == null || species2 == null || species3 == null) throw new NullPointerException("Cannot set null Species in Potential3");
		if(!(parentPotential() instanceof PotentialMaster)) throw new RuntimeException("Error: Can set species only for potentials that apply at the molecule level.  Potential must have PotentialMaster as parent");
		((PotentialMaster)parentPotential()).setSpecies(this, species);
    }
    /**
     * Returns an array of length 2 with the species to which this potential applies.
     * Returns null if no species has been set, which is the case if the potential
     * is not describing interactions between molecule-level Atoms.
     */
    public Species[] getSpecies() {
        if(species1 == null) return null;
        else return new Species[] {species1, species2, species3};
    }

    public void setIterator(Atom3Iterator iterator) {
        this.iterator = iterator;
    }
    public Atom3Iterator iterator() {return iterator;}
    
}//end of Potential3



