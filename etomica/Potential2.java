package etomica; 

/**
 * Potential acting on a pair of atoms or atom groups.  The Potential2Group
 * subclass of this potential is used to describe the inter-group interactions.
 * The Potential2Group would contain other Potential2 (or higher-body) instances
 * that describe the specific interactions between the atoms of the group.
 *
 * @author David Kofke
 */
public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:01.07.03/"+Potential.VERSION;
    
    protected AtomPairIterator iterator;
    protected AtomPairIterator iterator1;
    protected AtomPairIterator iteratorA;
    private Species species1, species2;
    
    public final PotentialTruncation truncation;

    public Potential2(PotentialGroup parent) {
        super(parent);
        iterator = new ApiGeneral(parentSimulation().space());
        iterator1 = new ApiIntergroup1A(parentSimulation());
        iteratorA = new ApiIntergroupAA(parentSimulation());
        if(Default.TRUNCATE_POTENTIALS) {//can't use other constructor because of "this" in constructor of PotentialTruncationSimple
            truncation = new PotentialTruncationSimple(this, Default.POTENTIAL_CUTOFF_FACTOR * Default.ATOM_SIZE);
        } else {
            truncation = PotentialTruncation.NULL;
        }
    }
    public Potential2(PotentialGroup parent, PotentialTruncation trunc) {
        super(parent);
        iterator = new ApiGeneral(parentSimulation().space());
        truncation = trunc;
    }
    
    public final void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential2Calculation) ) return;
        iterator = (id.atomCount() == 0) ? iteratorA : iterator1;
        iterator.reset(id);
        calculate2(id, (Potential2Calculation)pc);
    }//end of calculate

    public abstract void calculate2(IteratorDirective id, Potential2Calculation pc);
    public abstract double energy(AtomPair pair);

    public Potential set(Atom a) {return set(a,a);}
    //need to find a better way to set basis 
    //problem is that iterator is selected in calculate method, after set is called
    public Potential set(Atom a1, Atom a2) {
  //      iterator.setBasis(a1, a2); 
        iterator1.setBasis(a1, a2); 
        iteratorA.setBasis(a1, a2); 
        return this;}
    
    public Potential set(Atom[] atoms) {
        switch (atoms.length) {
            case 1: return set(atoms[0], atoms[0]);
            case 2: return set(atoms[0], atoms[1]);
            default: throw new IllegalArgumentException("Wrong number of atoms in Potentia2.set");
        }
    }
    public Potential set(SpeciesMaster s) {
        if(species1 == null || species2 == null) 
            throw new IllegalStateException("Attempt to set SpeciesMaster in Potential2 before species are specified");
  //      if(species1 != null) {//if species were previously set, use them as basis
 //           iterator.setBasis(species1.getAgent(s), species2.getAgent(s));
            iterator1.setBasis(species1.getAgent(s), species2.getAgent(s));
            iteratorA.setBasis(species1.getAgent(s), species2.getAgent(s));
   /*     } else {
            iterator.setBasis(s,s); //otherwise use speciesMaster as basis
            iterator1.setBasis(s,s); //otherwise use speciesMaster as basis
            iteratorA.setBasis(s,s); //otherwise use speciesMaster as basis
        }*/
        return this;
      //  iterator.setBasis((species != null) ? species.getAgent(s.parentPhase()) : s);
    }
    
    public void setSpecies(Species s1, Species s2) {
        if(s1 == null || s2 == null) throw new NullPointerException("Cannot set null Species in Potential2");
        species1 = s1;
        species2 = s2;
        if(species1 == species2) {
            iterator1 = new ApiIntraspecies1A(parentSimulation());
            iteratorA = new ApiIntraspeciesAA(parentSimulation());
        } else {
            iterator1 = new ApiInterspecies1A(parentSimulation());
            iteratorA = new ApiInterspeciesAA(parentSimulation());
        }
//        if(speciesMaster != null) set(speciesMaster);
    }

    public void setIterator(AtomPairIterator iterator) {
        this.iterator = iterator;
        iterator1 = iterator;
        iteratorA = iterator;
    }
    public AtomPairIterator iterator() {return iterator;}
    
    public void setSpecies(Species[] species) {
        switch (species.length) {
            case 1: setSpecies(species[0], species[0]);
                    break;
            case 2: setSpecies(species[0], species[1]);
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
        if(species1 == null) return null;
        else return new Species[] {species1, species2};
    }
    
    /**
     * Accessor method for potential cutoff implementation.
     */
    public PotentialTruncation getTruncation() {return truncation;}
    
}//end of Potential2



