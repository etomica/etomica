package etomica;

/**
 * Collection of zero-body (Potential0) potentials.  A set(Phase) method call
 * is passed to the contained potentials, and a calculate(IteratorDirective, PotentialCalculation)
 * method call is passed on to these potentials, which then perform the given
 * calculation using the previously-specified phase.
 *
 * @author David Kofke
 */

public class Potential0Group extends Potential0 implements PotentialGroup {
    
    private final PotentialCalculationEnergySum energy = null;//new PotentialCalculation.EnergySum();
    protected PotentialLinker first;
    protected Phase phase;
    
    public Potential0Group() {
        this(Simulation.instance.hamiltonian.potential);
    }
    public Potential0Group(PotentialGroup parent) {
        super(parent);
    }
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if(phase == null) return;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
            link.potential.calculate(id, pc);
        }//end for
    }//end calculate
    
    /**
     * Identifies the phase that is used for any subsequent property calculations
     * (until next call to this method).  Passes on the specification to any potential0
     * instances that have been added to this group.  If given phase is the same
     * as the one most recently specified, no action is taken by method.
     */
    public Potential set(Phase phase) {
        if(this.phase == phase) return this;
        this.phase = phase;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            ((Potential0)link.potential).set(phase);
        }//end for
        return this;
    }

    /**
     * Adds potential to group and sets its phase to the one most recently
     * passed to the set(Phase) method.  If given potential is null, no
     * action is performed.
     */
    public void addPotential(Potential potential) {
        if(potential == null) return;
        if( !(potential instanceof Potential0) ) throw new IllegalArgumentException("Error:  Attempt to add to Potential0Group a potential that is not an instance of Potential0");
        first = new PotentialLinker(potential, first);
        if(phase != null) ((Potential0)potential).set(phase);
    }

    /**
     * Not implemented.
     */
    public double energy() {
        throw new RuntimeException("energy method not implemented in Potential0Group");
    }

}//end Potential0Group
    