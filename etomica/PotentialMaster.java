package etomica;

/**
 * Master potential that oversees all other potentials in the Hamiltonian.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation 
 * on to the contained potentials.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 8/13/02 (DAK) added removePotential method
  */
public final class PotentialMaster implements PotentialGroup, java.io.Serializable {
    
    public String getVersion() {return "PotentialMaster:02.08.13";}

    private SpeciesMaster speciesMaster;
    private Potential0GroupLrc lrcMaster;
    private final Simulation parentSimulation;
    private PotentialLinker first, last;

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
        
    /**
     * Adds the potential to the group.  Normally invoked in the constructor of the
     * given potential.
     */
    public void addPotential(Potential potential) {
        //Set up to evaluate zero-body potentials last, since they may need other potentials
        //to be configured for calculation (i.e., iterators set up) first
        if(potential instanceof Potential0 && last != null) {//put zero-body potential at end of list
            last.next = new PotentialLinker(potential, null);
            last = last.next;
        } else {//put other potentials at beginning of list
            first = new PotentialLinker(potential, first);
            if(last == null) last = first;
        }
        //if phase was set previously, subsequent call to set to same phase will be
        //ignored (and not passed on to new potential).  Seting speciesMaster to null
        //here will prevent ignoring of new call to set phase 
        speciesMaster = null;
    //    if(speciesMaster != null) {
    //        System.out.println("Warning: adding potential after phase was set for PotentialMaster");
    //        System.out.println("May lead to error");
    //        speciesMaster = null;
    //    }
//        if(speciesMaster != null) potential.set(speciesMaster); can't do this because potential is still executing its constructor
    }
    
    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public void removePotential(Potential potential) {
        PotentialLinker previous = null;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential == potential) {//found it
                if(previous == null) first = link.next;  //it's the first one
                else previous.next = link.next;          //it's not the first one
                if(link == last) last = previous; //removing last; this works also if last was also first (then removing only, and set last to null)
                return;
            }//end if
            previous = link;
        }//end for
    }//end removePotential

    /**
     * Sets the basis for iteration of atoms by all potentials.
     * The given parameter must be an instance of SpeciesMaster.
     * No action is taken if the given species master is the same
     * as the one given in the previous call (unless another
     * potential was added in the interim).
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
    