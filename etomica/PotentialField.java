package etomica; 

/**
 * Basic methods needed to describe the action of an external field (a single-body potential) on an atom.
 * A field is introduced by adding it to a Phase.  See PotentialFieldGravity for a main method demonstrating
 * the use of a field.
 */
public abstract class PotentialField extends PotentialAbstract implements java.io.Serializable {
  
    private PotentialField next;
    protected Phase phase;
    private Atom.Iterator atomIterator;
    protected Maker maker;
    
    public PotentialField(Phase p) {
        this(p, p.iteratorFactory().makeAtomIterator()); //default atom iterator is all atoms
    }
    public PotentialField(Phase p, Atom.Iterator iter) {
        phase = p;
        atomIterator = iter; 
    }
    
    /**
     * Object that made this field.
     */
    public Maker maker() {return maker;}
    
    public final PotentialField nextField() {return next;}
    public final void setNextField(PotentialField field) {next = field;}
    
    public final Phase phase() {return phase;}
    
    /**
     * Iterator for atoms under the influence of this field.
     * Default is iterator of all atoms in the phase
     */
    public Atom.Iterator getAffectedAtoms() {return atomIterator;}
    public void setAffectedAtoms(Atom.Iterator iterator) {atomIterator = iterator;}
        
    /**
     * Returns the energy due to the interaction of the atom with the field.
     */
    public abstract double energy(Atom atom);
    
    /**
     * Returns the total energy of the field with all affected atoms in the phase
     */
    public double energy() {
        double sum = 0.0;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            sum += energy(atomIterator.next());
        }
        return sum;
    }
          
    
    //***************** end of methods for PotentialField class *****************//
    
    /**
     * Interface for a class that can make a PotentialField object.
     */
    public interface Maker {
        public PotentialField makePotentialField(Phase p);
    }
    
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



