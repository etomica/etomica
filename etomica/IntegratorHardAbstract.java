package etomica;

/**
 * Parent class of all integrators for hard potentials.
 * Integrates equations of motion through time by advancing atoms from one collision to the
 * next.  Management of sequence of collisions is delegated to the subclass.
 * Determination of time of collision and implementation of collision
 * dynamics is handled by the potential between the atoms, which therefore must
 * implement PotentialHard.
 *
 * @see PotentialHard
 * @author David Kofke
 *
 */
 
public abstract class IntegratorHardAbstract extends IntegratorMD {

    public static String VERSION = "IntegratorHardAbstract:01.06.14/"+IntegratorMD.VERSION;

    //convenience handle to the agent holding information about the next collision
    protected Agent nextCollider;
    //iterators for looping through atoms
    protected Atom.Iterator upAtomIterator;
    //first of a linked list of objects (typically meters) that are called each time a collision is processed
    protected CollisionListenerLinker collisionListenerHead = null;
    //time elapsed since reaching last timestep increment
    private double timeIncrement = 0.0;
    protected AtomPair atomPair;
    protected Potential.Hard potential;
                
    public IntegratorHardAbstract(Simulation sim) {
        super(sim);
    }

	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
        upAtomIterator = factory.makeAtomIteratorUp();
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
            System.out.println("NextCollider: "+nextCollider.atom);
            System.out.println("Collision partner: "+nextCollider.collisionPartner);
            System.out.println("Collision potential: "+nextCollider.collisionPotential);
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
        double collisionTimeStep = nextCollider.collisionTime();
        if(tStep < collisionTimeStep) {
            advanceAcrossTimeStep(tStep);
            timeIncrement = 0.0;
            if(isothermal) {
                scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
            }
        }
        else {
            double tStepNew = tStep - collisionTimeStep;
            advanceAcrossTimeStep(collisionTimeStep);
            processCollision();
            for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
                cll.listener.collisionAction(atomPair, potential); //atompair is not updated?
            }
            timeIncrement += collisionTimeStep;
            doStep(tStepNew);
        }
    }//end of doStep


    /**
    * Identifies atom that collides next, given that all collision times are known.
    */
    protected abstract void findNextCollider();

    /**
    * Applies collision dynamics to colliders and updates collision time/partners.
    */
    protected abstract void processCollision();

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
        upAtomIterator.reset();
        while(upAtomIterator.hasNext()) {
            Atom a = upAtomIterator.next();
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
     * Processes information about collisions between atoms.  A collision handler is
     * passed to the Potential and is informed about each collision detected by the
     * potential.  What the handler does with this information depends on how the integrator
     * is designed to manage the sequence of collisions.
     */
    public abstract static class CollisionHandler {
        public Potential.Hard potential;
        public final void setPotential(Potential.Hard p) {potential = p;}
        public abstract void setAtom(Atom a);
        public abstract void addCollision(Atom atom2, double collisionTime);
    }//end of CollisionHandler
    
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
        public Atom atom;
 //       public double time0;  //time since last collision
        private double collisionTime = Double.MAX_VALUE; //time to next collision
        Atom collisionPartner;  //next atom scheduled for collision by atom containing this Agent
        Potential.Hard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        PotentialField.Hard collisionPotentialField;
        boolean isFieldCollision;
        
        public Agent(Atom a) {atom = a;}
        
        public void resetCollision() {collisionTime = Double.MAX_VALUE;}
        
    /**
     * Sets parameters associated with next two-body collision of this atom with another atom.
     *
     * @param time    time to collision of this agent's atom with an atom uplist of it
     * @param partner the atom this one will collide with next
     * @param p       the potential for interactions between this atom and its collision partner
     */
        public final void setCollision(double time, Atom partner, Potential.Hard p) {
            isFieldCollision = false;
            collisionTime = time;
            collisionPartner = partner;
            collisionPotential = p;
        }

    /**
     * Sets parameters associated with next one-body collision of the atom with a hard field.
     *
     * @param time    time to collision of this agent's atom with an atom uplist of it
     * @param p       the potential for interactions between this atom and the field
     */
        public final void setCollision(double time, PotentialField.Hard p) {
            isFieldCollision = true;
            collisionTime = time;
            collisionPotentialField = p;
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
            if(a == nextCollider.atom) a.setColor(colliderColor);
            else if(a == nextCollider.collisionPartner) a.setColor(partnerColor);
            else a.setColor(baseColor);
        }
    }//end of HighlightColliders

    public interface CollisionListener {
        public void collisionAction(AtomPair pair, Potential.Hard p);
    }
    
}//end of IntegratorHardAbstract

