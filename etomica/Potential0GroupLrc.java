package etomica;

/**
 * Group that contains potentials used for long-range correction.
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
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if(phase == null || !phase.isLrcEnabled() || !id.includeLrc) return;
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(id.excludes(link.potential)) continue; //see if potential is ok with iterator directive
            link.potential.calculate(id, pc);
        }//end for
    }//end calculate
    
}