package etomica;

/**
 * Master potential that oversees all other potentials in the Hamiltonian.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation 
 * on to the contained potentials.
 *
 * @author David Kofke
 */
public final class PotentialMaster implements PotentialGroup, java.io.Serializable {
    
    public String getVersion() {return "PotentialMaster:01.07.23";}

    private SpeciesMaster speciesMaster;
    private Potential0GroupLrc lrcMaster;
    private final Simulation parentSimulation;
    private PotentialLinker first;

    public PotentialMaster(Simulation sim) {
        parentSimulation = sim;
    }
    
    public Simulation parentSimulation() {return parentSimulation;}
    
    //should build on this to do more filtering of potentials based on directive
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
            link.potential.calculate(id, pc);
        }//end for
    }//end calculate
        
    public final PotentialCalculation.Sum calculate(IteratorDirective id, PotentialCalculation.Sum pa) {
        this.calculate(id, (PotentialCalculation)pa);
        return pa;
    }
        
    //this method is called in the constructor of the given potential
    /**
     * Adds the potential to the group 
     */
     /*and calls its set(SpeciesMaster)
     * method if this.set(SpeciesMaster) was previously called with a
     * non-null argument.
     */
    public void addPotential(Potential potential) {
        first = new PotentialLinker(potential, first);
        if(speciesMaster != null) {
            System.out.println("Warning: adding potential after phase was set for PotentialMaster");
            System.out.println("May lead to error");
            speciesMaster = null;
        }
//        if(speciesMaster != null) potential.set(speciesMaster); can't do this because potential is still executing its constructor
    }

    /**
     * Sets the basis for iteration of atoms by all potentials.
     * The given parameter must be an instance of SpeciesMaster.
     * No action is taken if the given species master is the same
     * as the one given in the previous call.
     */
    public PotentialMaster set(SpeciesMaster sm) {
        if(sm == speciesMaster) return this;
        speciesMaster = sm;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            link.potential.set(speciesMaster);
        }//end for
        return this;
    }

    /**
     * Sets the potentials to iterate on atoms in the given phase.
     */
    public PotentialMaster set(Phase p) {
        return set(p.speciesMaster); 
    }
    
    /**
     * Returns the potential group that oversees the long-range
     * correction zero-body potentials.
     */
     public Potential0GroupLrc lrcMaster() {
        if(lrcMaster == null) lrcMaster = new Potential0GroupLrc(this);
        return lrcMaster;
     }

}//end of PotentialMaster
    