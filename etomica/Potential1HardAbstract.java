package etomica;

/**
 * Default abstract implementation of a one-body hard potential.
 */
public abstract class Potential1HardAbstract extends Potential1 implements Potential1.Hard {
    
    public static String VERSION = "Potential1HardAbstract:01.06.24/" + Potential1.VERSION;
    
    public Potential1HardAbstract(Simulation sim) {
        super(sim);
    }
    
    /**
    * Implements the collision dynamics.
    * The given atom is assumed to be at the point of collision.  This method is called
    * to change its momentum according to the action of the collision.  Extensions can be defined to
    * instead implement other, perhaps unphysical changes.
    */
    public abstract void bump(Atom atom);
    /**
    * Computes the time of collision of the given atom with the potential, 
    * assuming no intervening collisions.
    * Usually assumes free-flight between collisions.
    */ 
    public abstract double collisionTime(Atom atom);
    
    //abstract method energy(Atom) is inherited from Potential1 and must be
    //be defined in subclass

    public PotentialAgent makeAgent(Phase p) {
        return new Agent(this, p);
    }
        
    //Potential1HardAbstract.Agent
    public class Agent extends Potential1.Agent implements PotentialAgent.Hard {
        public Agent(PotentialAbstract potential, Phase phase) {
            super(potential, phase);
        }
        public void findCollisions(IteratorDirective id, 
                                    final IntegratorHardAbstract.CollisionHandler collisionHandler) {
            //collisions with one-body hard potentials are considered "uplist" of
            //all atoms, so no collisions with potential arise when doing a downlist iteration
            if(id.direction() == IteratorDirective.DOWN) return;
            
            collisionHandler.setPotential(this);
            iterator.reset(id);
            while(iterator.hasNext()) {
                Atom atom = iterator.next();
                double time = collisionTime(atom);
               //if(time >= 0) etc.  (maybe this would be more efficient)
                if(time < Double.MAX_VALUE) collisionHandler.addCollision(atom, time);
            }//end while
        }//end of findCollisions
    
        public void bump(IntegratorHardAbstract.Agent agent) {
            Potential1HardAbstract.this.bump(agent.atom);
        }
        
    }//end of Potential1HardAbstract.Agent
    
}//end of Potential1HardAbstract
    
    