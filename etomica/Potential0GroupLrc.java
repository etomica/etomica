package etomica;

/**
 * Group that contains potentials used for long-range correction.
 * One instance of this class is added to the PotentialManager of
 * a Simulation when the first Potential0Lrc class is instantiated.
 * All subsequently created Potential0Lrc classes are added to this group.
 *
 * @see Potential0Lrc
 * @author David Kofke
 */
 
public class Potential0GroupLrc extends Potential0Group {

    
/*    public Potential0GroupLrc() {
        this(Simulation.instance.hamiltonian.potential);
    }*/
    
    public Potential0GroupLrc(PotentialGroup parent) {
        super(parent);
    }
    
    /**
     * Performs given PotentialCalculation on all LRC potentials added to this group.
     * Checks that phase has been previously set, that it has lrcEnabled, and that
     * the given IteratorDirective has includeLrc set to true; if all are so, calculation
     * is performed.
     */
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if(phase == null || !phase.isLrcEnabled() || !id.includeLrc) return;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
            link.potential.calculate(id, pc);
        }//end for
    }//end calculate
    
}