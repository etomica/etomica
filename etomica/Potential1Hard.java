package etomica;
    
/**
 * Methods needed to describe the behavior of a hard one-body potential.  
 * A hard potential describes impulsive interactions, in which the energy undergoes a step
 * change at some point in the space.
 */
 
 /* History of changes
  * 08/26/02 (DAK) Modified calculate to handle cases of atomCount 0 or 1 in iterator directive
  */
public abstract class Potential1Hard extends Potential1Group implements PotentialHard {

    public Potential1Hard(PotentialGroup parent) {
        super(parent);
    }
   public void bump(AtomPair pair) {bump(pair.atom1());}
   
 /**
  * Implements the collision dynamics.
  * The given atom is assumed to be at the point of collision.  This method is called
  * to change its momentum according to the action of the collision.  Extensions can be defined to
  * instead implement other, perhaps unphysical changes.
  */
    public abstract void bump(Atom atom);

 /**
  * Computes the time of collision of the given atom with the hard potential, assuming no intervening collisions.
  * Usually assumes free-flight between collisions
  */ 
    public abstract double collisionTime(Atom atom);
            
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential1Calculation) ) return;
        if(id.atomCount() == 0) {
            iterator.reset(id);
            ((Potential1Calculation)pc).calculate(iterator, this); 
        } else {
            singletIterator.setAtom(id.atom1());
            singletIterator.reset();
            ((Potential1Calculation)pc).calculate(singletIterator, this); 
        }            
    }//end of calculate
    
    private AtomIteratorSinglet singletIterator = new AtomIteratorSinglet();

}  //end of Potential1Hard

