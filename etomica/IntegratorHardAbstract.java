package etomica;

/**
 * Parent class of all integrators for hard potentials.
 * Integrates equations of motion through time by advancing atoms from one collision to the
 * next.  Management of sequence of collisions is delegated to the subclass.
 * Determination of time of collision and implementation of collision
 * dynamics is handled by the potential between the atoms, which therefore must
 * have an Agent that implements PotentialAgent.Hard.
 * <br>
 * This class is not thread safe.
 *
 * @author David Kofke
 *
 */
 
 //thread safety of this class is compromised by the
 //Space.Vector methods used in the periodCollisionTime method of Agent
 
public abstract class IntegratorHardAbstract extends IntegratorMD {

    public static String VERSION = "IntegratorHardAbstract:01.07.03/"+IntegratorMD.VERSION;

    //handle to the integrator agent holding information about the next collision
    protected IntegratorHardAbstract.Agent colliderAgent;
    //iterators for looping through atoms
    protected AtomIterator atomIterator;
    //first of a linked list of objects (typically meters) that are called each time a collision is processed
    protected CollisionListenerLinker collisionListenerHead = null;
    //time elapsed since reaching last timestep increment
    private double timeIncrement = 0.0;
    protected PotentialAgent phasePotential;
    private AtomPair atomPair;
    Space.Vector c3;
                
    public IntegratorHardAbstract(Simulation sim) {
        super(sim);
        c3 = sim.space().makeVector();
        c3.E(3.0);
        
    }//end of constructor

	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
	    super.makeIterators(factory);
        atomIterator = factory.makeAtomIterator();
    }
    
          //need to modify to handle multiple-phase issues
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        phasePotential = p.potential();
        atomPair = new AtomPair(p);
        return true;
    }

    /** 
     * Steps all atoms across time interval timeStep, handling all intervening collisions.
     */
    public void doStep() {
        try {doStep(timeStep);}
        catch(StackOverflowError e) {
            System.out.println();
            System.out.println("Stack overflow error encountered when attempting to complete time step");
            System.out.println("This occurs when the simulated system becomes jammed, and is unable to clear all collisions and move to the next time step.");
            System.out.println("nextCollider: "+colliderAgent.atom());
            System.out.println("Collision partner: "+colliderAgent.collisionPartner());
            System.out.println("Collision potential: "+colliderAgent.collisionPotential);
            System.exit(1);
        }
    }//end of doStep
    
    /**
    * Steps all atoms across time interval tStep, handling all intervening collisions.
    * This method is called recursively, each time with a timestep equal to the amount of 
    * time between the latest collision and the endpoint of the timestep identified in its
    * highest-level call
    */
    public void doStep(double tStep) {
        double collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.MAX_VALUE;
        if(tStep < collisionTimeStep) {
            advanceAcrossTimeStep(tStep);
            timeIncrement = 0.0;
            if(isothermal) {
                scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
            }
        }
        else {
            double tStepNew = tStep - collisionTimeStep;
 //           System.out.println(colliderAgent.atom().toString() +" "+ colliderAgent.collisionPartner().toString());
            advanceAcrossTimeStep(collisionTimeStep);//if needing more flexibility, make this a separate method-- advanceToCollision(collisionTimeStep)
            atomPair.reset(colliderAgent.atom(), colliderAgent.collisionPartner());
            
            colliderAgent.collisionPotential.bump(atomPair);
            
            for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
                cll.listener.collisionAction(colliderAgent);
            }
            updateCollisions();
            findNextCollider(); //this sets colliderAgent for the next collision
            timeIncrement += collisionTimeStep;
            doStep(tStepNew);
        }
    }//end of doStep


    /**
    * Identifies atom that collides next, given that all collision times are known.
    */
    protected abstract void findNextCollider();

    /**
    * Updates collision times/partners.
    */
    protected abstract void updateCollisions();

    /**
    * Advances all atom coordinates by tStep, without any intervening collisions.
    */
    protected abstract void advanceAcrossTimeStep(double tStep);

    /**
    * The simulation time elapsed since the start of the integration.
    * Cannot be reset to zero.  
    * Modified from superclass method to return correct value if in the middle
    * of processing collisions between time steps.
    */
    public double elapsedTime() {return super.elapsedTime() + timeIncrement;}

              
    /**
    * Crude method to enforce constant-temperature constraint
    * Scales momenta of all atoms by a constant factor so that 
    * phase adheres to setpoint temperature.  Updates collision times appropriately.
    */
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            a.momentum().TE(s); //scale momentum
            ((Agent)a.ia).collisionTime *= rs;
        }
    }
        
    /**
     * Registers an object that implements the CollisionListener interface.
     * This causes the collisionAction method of the object to be called after each collision.
     * If the listener is already in the list, no action is taken to add it again.
     */
    public void addCollisionListener(CollisionListener cl) {
        //make sure listener is not already in list
        for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
            if(cll.listener == cl) return;
        }
        //OK, not in list, now register it by adding it to the beginning of the list
        collisionListenerHead = new CollisionListenerLinker(cl, collisionListenerHead);
    }
    /**
     * De-registers an object from the list of collision listeners.
     * No exception is generated if the listener is not in the list; 
     * method simply returns without taking any action.
     */
    public void removeCollisionListener(CollisionListener cl) {
        if(collisionListenerHead == null) return;   //list is empty; do nothing
        if(cl == collisionListenerHead.listener) {  //listener is first in list
            collisionListenerHead = collisionListenerHead.next;
        }
        else {                                     //not first, so look for listener later in list
            CollisionListenerLinker previous = collisionListenerHead;
            for(CollisionListenerLinker cll=collisionListenerHead.next; cll!=null; cll=cll.next) {
                if(cll.listener == cl) {           //found it; 
                    previous.next = cll.next;      //remove it from list
                    break;
                }
            }
        }
    }//end of removeCollisionListener
    
    /**
     * Class used to construct a linked list of collision listeners
     */
    private class CollisionListenerLinker {
        CollisionListenerLinker next;
        final CollisionListener listener;
        CollisionListenerLinker(CollisionListener listener, CollisionListenerLinker next) {
            this.listener = listener;
            this.next = next;
        }
    }//end of CollisionListenerLinker

 /**
  * Agent defined by this integrator.
  * Holds information about the time to next collision (considering only atoms
  * uplist of the atom), the collisionPartner (atom that this one will be colliding with),
  * and a handle to the potential governing the interactions of this atom and its collision partner
  * (this is kept for convenience, so it won't have to be determined again when the collision
  * is processed).
  */
  //Do not use encapsulation since the fields are for referencing by the integrator
    public static class Agent implements Integrator.Agent {  //need public so to use with instanceof
        public Atom atom, collisionPartner;
        public double collisionTime = Double.MAX_VALUE; //time to next collision
        public PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        public boolean movedInList = false;
        
        public Agent(Atom a) {atom = a;}
        
        public final Atom atom() {return atom;}
        public final Atom collisionPartner() {return collisionPartner;}
        
        public void resetCollision() {
            collisionTime = periodCollisionTime();
            collisionPotential = Potential2Hard.NULL;
            collisionPartner = null;
        }
        
        //time to "collision" to update colliders for periodic boundaries
        private Space.Vector c3;
        protected double periodCollisionTime() {
            Space.Boundary boundary = atom.parentPhase().boundary();
            if(boundary instanceof Space.Boundary.Periodic) {
                if(!(atom.type instanceof AtomType.Disk)) {return Double.MAX_VALUE;}
                Space.Vector p = atom.coordinate().momentum();
                Space.Vector dim = boundary.dimensions();
                double diameter = ((AtomType.Disk)atom.type).diameter();
                return 0.5*atom.mass()*dim.D(p).abs().min(); //0.5*m*min of (dim.x/p.x, dim.y/p.y, etc.)
          //      System.out.println(xx);
          //     return xx;
          //      return 0.5*atom.mass()*dim.D(p).abs().min(); //0.5*m*min of (dim.x/p.x, dim.y/p.y, etc.)
                //assumes range of potential is .le. diameter, simulation box is square (or x is smaller dimension)
            //    return 0.5*(dimensions.y-1.0001*diameter)/(a.rm()*Math.sqrt(p.squared()));
            }
            else return Double.MAX_VALUE;
        }
        
    /**
     * Sets parameters associated with next two-body collision of this atom with another atom.
     *
     * @param time    time to collision of this agent's atom with an atom uplist of it
     * @param partner the atom this one will collide with next
     * @param p       the potential for interactions between this atom and its collision partner
     */
        public final void setCollision(double time, Atom partner, PotentialHard p) {
            collisionTime = time;
            collisionPartner = partner;
            collisionPotential = p;
        }


    /**
     * Decreases the recorded time to collision of this atom
     * This action is performed when the atom is advanced without a collision
     */
        public final void decrementCollisionTime(double interval) {
            collisionTime -= interval;
//            time0 += interval;
        }
    /**
     * Accessor method for the time to next collision of this atom
     */
        public final double collisionTime() {return collisionTime;}
    }
    
    /**
     * This colorScheme acts to color differently the two atoms that are scheduled to collide next.
     * Highlight colors are specified by the colliderColor and partnerColor fields; all other
     * atoms are colored with the baseColor
     */
    public class HighlightColliders extends ColorScheme {
        
        public HighlightColliders() {super();}
        /**
         * Color applied to the downList atom of the colliding pair
         */
        public java.awt.Color colliderColor = java.awt.Color.red;
        /**
         * Color applied to the upList atom of the colliding pair
         */
        public java.awt.Color partnerColor = java.awt.Color.blue;
        /**
         * Applies the special colors to the colliding pair while coloring all other atoms with baseColor.
         */ 
        public void colorAtom(Atom a) {
            if(colliderAgent == null) a.setColor(baseColor);
            else if(a == colliderAgent.atom) a.setColor(colliderColor);
            else if(a == colliderAgent.collisionPartner) a.setColor(partnerColor);
            else a.setColor(baseColor);
        }
    }//end of HighlightColliders

    public interface CollisionListener {
        public void collisionAction(Agent colliderAgent);
    }
    
}//end of IntegratorHardAbstract

