package etomica; 

/**
 * Potential acting on a single atom group. This could be an external field acting
 * on a single atom, or could describe the interactions among the atoms in a group.
 */
public abstract class Potential1 extends PotentialAbstract {
  
    public static String VERSION = "Potential1:01.06.12/"+PotentialAbstract.VERSION;

    private Atom1Iterator iterator;
    
    public Potential1(Simulation sim) {
        super(sim);
    }
    
    public Atom1Iterator iterator() {return iterator;}
    
    /**
     * Returns the energy of the given atom group.
     */
    public abstract double energy(AtomGroup atom);
    
    /**
     * Returns the total energy of the field with all affected atoms in the phase
     */
    public double energy() {
        double sum = 0.0;
        iterator.reset();
        while(iterator.hasNext()) {
            sum += energy(iterator.next());
        }
        return sum;
    }
          
    
    //***************** end of methods for Potential1 class *****************//
    
    
    /**
    * Methods needed to describe the behavior of a hard field potential.  
    * A hard potential describes impulsive interactions, in which the energy undergoes a step
    * change at some point in the space.
    *
    * @see PotentialField.Soft
    */
         
    public interface Hard {
            
    /**
     * Returns the energy due to the interaction of the atom with the field.
     */
        public double energy(Atom atom);
    /**
    * Implements the collision dynamics.
    * The given atom is assumed to be at the point of collision.  This method is called
    * to change its momentum according to the action of the collision.  Extensions can be defined to
    * instead implement other, perhaps unphysical changes.
    */
        public void bump(Atom atom);
    /**
    * Computes the time of collision of the given atom with the external field, assuming no intervening collisions.
    * Usually assumes free-flight between collisions
    */ 
        public double collisionTime(Atom atom);
            
    }  //end of PotentialField.Hard

    /**
    * Methods needed to describe the behavior of a soft potential.  
    * A soft potential describes non-impulsive interactions, in which the energy at all points
    * has smooth, analytic behavior with no discontinuities.  
    *
    * @see PotentialField.Hard
    */
    public interface Soft {
        
        /**
        * Returns the energy due to the interaction of the atom with the field.
        */
        public double energy(Atom atom);

        /**
        * Force exerted by the field on the atom.
        *
        * @return the vector force exerted on the atom
        */
        public Space.Vector force(Atom atom);

    } //end of PotentialField.Soft
}



