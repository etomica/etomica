package etomica; 

/**
 * Potential acting on a pair of atom groups.
 */
public abstract class Potential2 extends PotentialAbstract {
  
    public static String VERSION = "Potential2:01.06.12/"+PotentialAbstract.VERSION;
    
    private AtomPair.Iterator iterator;
    
    public Potential2(Simulation sim) {
        super(sim);
    }
    
    public AtomPair.Iterator iterator() {return iterator;}
    
    /**
     * Returns the energy of the given atom group.
     */
    public abstract double energy(AtomPair pair);
    
/*    public double energy(AtomGroup group) {
        iterator.reset(group);
        return energySum();
    }
*/    
    public double energy(Atom atom) {
        iterator.reset(atom);
        return energySum();
    }
    
    /**
     * Returns the total energy of all affected atoms.
     */
    public double energy() {
        iterator.reset();
        return energySum();
    }
    
    private final double energySum() {
        double sum = 0.0;
        while(iterator.hasNext()) {
            sum += energy(iterator.next());
        }
        return sum;
    }
          
    
    //***************** end of methods for Potential2 class *****************//
    
    
    /**
    * Methods needed to describe the behavior of a hard potential.  
    * A hard potential describes impulsive interactions, in which the energy undergoes a step
    * change at some point in the space.
    *
    * @see PotentialField.Soft
    */
         
    public interface Hard {
            
    /**
     * Returns the energy due to the interaction of the atoms with each other
     */
        public double energy(AtomPair pair);
    /**
    * Implements the collision dynamics.
    * The given atoms are assumed to be at the point of collision.  This method is called
    * to change their momentum according to the action of the collision.  Extensions can be defined to
    * instead implement other, perhaps unphysical changes.
    */
        public void bump(AtomPair pair);
    /**
    * Computes the time of collision of the given atoms , assuming no intervening collisions.
    * Usually assumes free-flight between collisions
    */ 
        public double collisionTime(AtomPair pair);
            
    }  //end of Potential2.Hard

    /**
    * Methods needed to describe the behavior of a soft potential.  
    * A soft potential describes non-impulsive interactions, in which the energy at all points
    * has smooth, analytic behavior with no discontinuities.  
    *
    * @see PotentialField.Hard
    */
    public interface Soft {
        
        /**
        * Returns the energy due to the interaction of the atom with each other.
        */
        public double energy(AtomPair pair);

        /**
        * Force exerted atoms exert on each other.
        *
        * @return the vector force exerted on the atom
        */
        public Space.Vector force(AtomPair pair);

    } //end of PotentialField.Soft
}



