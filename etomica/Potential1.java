package etomica; 

/**
 * Potential acting on a single atom.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends PotentialAbstract {
  
    public static String VERSION = "Potential1:01.06.27/"+PotentialAbstract.VERSION;
    
    public Potential1(Simulation sim) {
        super(sim);
    }
    
    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
          
    
    //***************** end of methods for Potential1 class *****************//
    
    //Potential1.Agent
    public class Agent extends PotentialAgent {
        
        protected AtomIterator iterator;
        /**
         * @param potential The parent potential making this agent
         * @param phase The phase in which this agent will be placed
         */
        public Agent(PotentialAbstract potential, Phase phase) {
            super(potential, phase);
        }
        
        protected void makeIterator() {
            iterator = AtomIterator.NULL;
        }
            
        public void setIterator(AtomIterator iterator) {
            this.iterator = iterator;
        }
        public AtomIterator iterator() {return iterator;}
    
        public final PotentialAbstract parentPotential() {return Potential1.this;}
        
    
       /**
        * Returns the total energy of the potential with all affected atoms.
        */
        public double energy(IteratorDirective id) {
            iterator.reset(id);
            double sum = 0.0;
            while(iterator.hasNext()) {
                sum += Potential1.this.energy(iterator.next());
            }
            return sum;
        }
    }//end of Agent    
    
    /**
    * Methods needed to describe the behavior of a hard field potential.  
    * A hard potential describes impulsive interactions, in which the energy undergoes a step
    * change at some point in the space.
    */
    public interface Hard {

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
            
    }  //end of Potential1.Hard

    /**
    * Methods needed to describe the behavior of a soft potential.  
    * A soft potential describes non-impulsive interactions, in which the energy at all points
    * has smooth, analytic behavior with no discontinuities.  
    *
    * @see PotentialField.Hard
    */
    public interface Soft {
        
        /**
        * Force exerted by the field on the atom.
        *
        * @return the vector force exerted on the atom
        */
        public Space.Vector force(Atom atom);

    } //end of PotentialField.Soft
}



