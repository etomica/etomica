package etomica;

/**
 * Valence force-field potential, used in modeling of semiconductors.  Consists
 * of two-body and three-body parts.
 */
 
public class PotentialVFF extends Potential {
    
    private final _2BodyPart potential2 = new _2BodyPart();
    private final _3BodyPart potential3 = new _3BodyPart();

    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        iterator2.reset(id);  //reset for iteration over pairs of atoms
        ((Potential2Calculation)pc).calculate(iterator2, potential2); 

        iterator3.reset(id);
        ((Potential3Calculation)pc).calculate(iterator3, potential3); 
    }//end calculate
        
    //Sets the basis for iteration
    public abstract Potential set(Atom a);
    public abstract Potential set(Atom a1, Atom a2);
    public abstract Potential set(SpeciesMaster s);    
    
        
    private static final class _2BodyPart extends Potential2Soft {

        //not implemented
        public double energy(AtomPair pair) {
            return 0.0;
        }
        
        //not implemented
        public double virial(AtomPair pair) {
            return 0.0;
        }
        
        //not implemented
        public double hyperVirial(AtomPair pair) {
            return 0.0;
        }
        
        //not implemented
        public Space.Vector gradient(AtomPair pair) {
            return null;
        }
        
        /**
        * Integral used to evaluate correction to truncation of potential.
        */
        public double integral(double rC) {
            return 0.0;
        }
    }//end of _2BodyPart
    
    private static final class _3BodyPart extends Potential3 {
        public double energy(Atom3 atom3) {
        }
    }//end of _3BodyPart
}