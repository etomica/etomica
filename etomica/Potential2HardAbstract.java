package etomica;

/**
 * Default abstract implementation of a two-body hard potential.
 */
public abstract class Potential2HardAbstract extends Potential2 implements Potential2.Hard {
    
    public static String VERSION = "Potential2HardAbstract:01.06.24/" + Potential2.VERSION;
    
    public Potential2HardAbstract(Simulation sim) {
        super(sim);
    }
    
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
    
    //abstract method energy(AtomPair) is inherited from Potential2 and must be
    //be defined in subclass

    public PotentialAgent makeAgent(Phase p) {
        return new Agent(p);
    }
        
    //P2HardAbstract.Agent
    public class Agent extends Potential2.Agent implements PotentialAgent.Hard {
        private final AtomPair atomPair;
        public Agent(Phase p) {
            super(p);
            atomPair = new AtomPair(p);
        }
        public void findCollisions(IteratorDirective id, 
                                    final IntegratorHardAbstract.CollisionHandler collisionHandler) {
            collisionHandler.setPotential(this);
            iterator.reset(id);
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                double time = collisionTime(pair);
                if(time < Double.MAX_VALUE) collisionHandler.addCollision(pair, time);
                //if(time >= 0) etc.  (maybe this would be more efficient)
            }//end while
        }//end of findCollisions
    
        public void bump(IntegratorHardAbstract.Agent agent) {
            atomPair.reset(agent.atom, agent.collisionPartner);
            Potential2HardAbstract.this.bump(atomPair);
        }
        
    }//end of P2HardAbstract.Agent
    
}//end of P2HardAbstract
    
    