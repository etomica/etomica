package etomica; 

/**
 * Potential acting on a trio of atoms or atom groups.
 *
 * @author David Kofke
 */
public abstract class Potential3 extends Potential {
  
    public static String VERSION = "Potential2:01.08.06/"+Potential.VERSION;
    
    protected Atom3Iterator iterator;
    
    public Potential3(PotentialGroup parent) {
        super(3, parent);
//        iterator = new Atom3Iterator(parentSimulation().space());
    }
    
    public abstract double energy(Atom3 atom3);
    
    public void setSpecies(Species s1, Species s2, Species s3) {
		setSpecies(new Species[] {s1, s2, s3});
    }

    public void setSpecies(Species[] species) {
    	throw new etomica.exception.MethodNotImplementedException();
//		super.setSpecies(species);
//        switch (species.length) {
//            case 1: 
//                    break;
//            case 2:
//            		break;
//            case 3: 
//                    break;
//            default: throw new IllegalArgumentException("Wrong number of species given in Potential3");
//        }
    }

	public void setIterator(AtomSetIterator iterator) {
		if(iterator instanceof Atom3Iterator) this.setIterator((Atom3Iterator)iterator);
		else throw new IllegalArgumentException("Inappropriate type of iterator set for potential");
	}
	public void setIterator(Atom3Iterator iterator) {
		this.iterator = iterator;
	}
	public AtomSetIterator getIterator() {return iterator;}
    
}//end of Potential3



