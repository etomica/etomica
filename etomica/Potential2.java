package etomica; 

/**
 * Potential acting on a pair of atoms.
 *
 * @author David Kofke
 */
public abstract class Potential2 extends PotentialAbstract {
  
    public static String VERSION = "Potential2:01.06.24/"+PotentialAbstract.VERSION;
    
    public Potential2(Simulation sim) {
        super(sim);
    }
    
    /**
     * Returns the energy of the given atom pair.
     */
    public abstract double energy(AtomPair pair);
    
            
    //***************** end of methods for Potential2 class *****************//
    
    //Potential2.Agent
    public class Agent extends PotentialAgent {
        
        protected AtomPairIterator iterator;
        /**
         * @param p The phase in which this agent will be placed
         */
        public Agent(Phase p) {
            super(p);
//            iterator = new AtomPairIterator(p);
        }
            
        public void setIterator(AtomPairIterator iterator) {
            this.iterator = iterator;
        }
        public AtomPairIterator iterator() {return iterator;}
    
        public final PotentialAbstract parentPotential() {return Potential2.this;}
        
        public double energy(IteratorDirective id) {
            iterator.reset(id);
            double sum = 0.0;
            while(iterator.hasNext()) {
                sum += Potential2.this.energy(iterator.next());
            }
            return sum;
        }
    }//end of Agent    

    //Potential2.Hard
    public interface Hard {
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
    }//end of Hard

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



