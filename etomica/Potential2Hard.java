package etomica;

/**
 * Parent class for all hard pair potentials.
 *
 * @author David Kofke
 */
public abstract class Potential2Hard extends Potential2 implements PotentialHard {
    
    public Potential2Hard(PotentialGroup parent) {
        super(parent);
    }
    public abstract double energy(AtomPair pair);
    /**
     * Implements the collision dynamics.
     * The given atoms are assumed to be at the point of collision.  This method is called
     * to change their momentum according to the action of the collision.  Extensions can be defined to
     * instead implement other, perhaps unphysical changes.
     */
    public abstract void bump(AtomPair pair);
    
    /**
     * Computes the time of collision of the given atoms , assuming no intervening collisions.
     * Usually assumes free-flight between collisions
     */ 
    public abstract double collisionTime(AtomPair pair);
    
    public final void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential2Calculation) ) return;
        iterator = (id.atomCount() == 0) ? iteratorA : iterator1;
        iterator.reset(id);

        ((Potential2Calculation)pc).calculate(iterator, this); 
   ///     iterator.allPairs(((PotentialCalculationEnergySum)pc).getAtomPairCalculation(this)); 
    }//end of calculate

}//end of Potential2Hard

    