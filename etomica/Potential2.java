package etomica; 

/**
 * Potential acting on a pair of atoms or atom groups.
 *
 * @author David Kofke
 */
public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:01.07.03/"+Potential.VERSION;
    
    protected AtomPairIterator iterator;
    private Species species1, species2;
    
    public Potential2(PotentialGroup parent) {
        super(parent);
        iterator = new AtomPairIterator(parentSimulation().space());
    }
    
    public abstract double energy(AtomPair pair);

    public Potential set(Atom a) {return set(a,a);}
    public Potential set(Atom a1, Atom a2) {iterator.setBasis(a1, a2); return this;}
             
    public Potential set(SpeciesMaster s) {
        if(species1 != null) {//if species were previously set, use them as basis
            iterator.setBasis(species1.getAgent(s), species2.getAgent(s));
        }
        else iterator.setBasis(s,s); //otherwise use speciesMaster as basis
        return this;
      //  iterator.setBasis((species != null) ? species.getAgent(s.parentPhase()) : s);
    }
    
    public void setSpecies(Species s1, Species s2) {//should throw exception if either are null
        species1 = s1;
        species2 = s2;
    }

    public void setIterator(AtomPairIterator iterator) {
        this.iterator = iterator;
    }
    public AtomPairIterator iterator() {return iterator;}
    
}//end of Potential2



