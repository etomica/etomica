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

    public Potential set(Atom a) {return set(a,a,a);}
    public Potential set(Atom a1, Atom a2) {return null;}//throw exception?
    public Potential set(Atom a1, Atom a2, Atom a3) {iterator.setBasis(a1, a2, a3); return this;}
             
    public Potential set(SpeciesMaster s) {
        if(species1 != null) {//if species were previously set, use them as basis
            iterator.setBasis(species1.getAgent(s), species2.getAgent(s), species3.getAgent(s));
        }
        else iterator.setBasis(s,s,s); //otherwise use speciesMaster as basis
        return this;
      //  iterator.setBasis((species != null) ? species.getAgent(s.parentPhase()) : s);
    }
    
    public void setSpecies(Species s1, Species s2, Species s3) {//throw exception if either are null
        species1 = s1;
        species2 = s2;
        species3 = s3;
    }

    public void setIterator(Atom3Iterator iterator) {
        this.iterator = iterator;
    }
    public Atom3Iterator iterator() {return iterator;}
    
}//end of Potential3



