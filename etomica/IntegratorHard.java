package etomica;
import java.util.Observable;
import java.awt.Color;
import etomica.units.*;

/**
 * General-purpose integrator for hard potentials.
 * Integrates equations of motion through time by advancing atoms from one collision to the
 * next.  Determination of time of collision and implementation of collision
 * dynamics is handled by the potential between the atoms, which therefore must
 * implement PotentialHard.  Each atom keeps in its Agent (obtained from this integrator) the shortest
 * time to collision with any of the atoms uplist of it (as defined by the iterator from 
 * the phase).
 *
 * @see PotentialHard
 * @author David Kofke
 *
 */
 
 // needs revision to recognize PotentialField classes
 
public class IntegratorHard extends IntegratorMD implements EtomicaElement {

//convenience handle to the agent holding information about the next collision
private Agent nextCollider;
//iterators for looping through atoms
private AtomPair.Iterator downPairIterator;
private AtomPair.Iterator upPairIterator;
private Atom.Iterator upAtomIterator;
private AtomPair atomPair;
private Potential.Hard spacePotential;
//first of a linked list of objects (typically meters) that are called each time a collision is processed
private CollisionListenerLinker collisionListenerHead = null;
//time elapsed since reaching last timestep increment
private double timeIncrement = 0.0;

//boolean bb = true;  used in debugging
            
public IntegratorHard() {
    this(Simulation.instance);
}
public IntegratorHard(Simulation sim) {
    super(sim);
}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Collision-based molecular dynamics simulation of hard potentials");
        return info;
    }

/**
 * @deprecated; use addPhase instead
 */
 public void registerPhase(Phase p) {addPhase(p);}
 
 /**
 * Identifies to this integrator a phase containing the atoms it is integrating.
 */
 public boolean addPhase(Phase p) {
    if(!super.addPhase(p)) return false;
    atomPair = new AtomPair(p);
//    spacePotential = (Potential.Hard)Simulation.space().makePotential(p);
    return true;
}

	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
        upPairIterator = factory.makeAtomPairIteratorUp();
        downPairIterator = factory.makeAtomPairIteratorDown();
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
}
/**
 * Steps all atoms across time interval tStep, handling all intervening collisions.
 * This method is called recursively, each time with a timestep equal to the amount of 
 * time between the latest collision and the endpoint of the timestep identified in its
 * highest-level call
 */
public void doStep(double tStep) {
    if(tStep < nextCollider.getCollisionTime()) {
        advanceAcrossTimeStep(tStep);
        timeIncrement = 0.0;
        if(isothermal) {
            scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
        }
//        debugMethod();
        return;
    }

    double collisionTimeStep = nextCollider.getCollisionTime();
    double tStepNew = tStep - collisionTimeStep;
    advanceToCollision();
    timeIncrement += collisionTimeStep;
    doStep(tStepNew);
    return;
}

//debugging tools
// Colors atoms according to whether they are in the up or down neighborlist of the first atom
    private void debugMethod() {
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            a.setColor(Color.black);
        }
        Atom central = firstPhase.firstAtom();
//        central.setColor(((Space2DCell.Coordinate)central.coordinate).cell.color);
        upPairIterator.reset(central);
        while(upPairIterator.hasNext()) {
            AtomPair pair = upPairIterator.next();
            pair.atom2.setColor(Color.blue);
        }
        downPairIterator.reset(central);
        while(downPairIterator.hasNext()) {
            AtomPair pair = downPairIterator.next();
            pair.atom2.setColor(Color.green);
        }
        /*
        
        upAtomIterator.reset(central);
        while(upAtomIterator.hasNext()) {
            Atom atom = upAtomIterator.next();
            atom.setColor(Color.blue);
        }
        downAtomIterator.reset(central);
        while(downAtomIterator.hasNext()) {
            Atom atom = downAtomIterator.next();
            atom.setColor(Color.green);
        }        
        */
//        central.setColor(((Space2DCell.Coordinate)central.coordinate).cell.color);
        central.setColor(Color.red);
    }


/**
 * Loops through all atoms to identify the one with the smallest value of collisionTime
 * Collision time is obtained from the value stored in the Integrator.Agent from each atom.
 */
protected void findNextCollider() {
    //find next collision pair by looking for minimum collisionTime
    double minCollisionTime = Double.MAX_VALUE;
    upAtomIterator.reset();
    while(upAtomIterator.hasNext()) {
        Agent ia = (Agent)upAtomIterator.next().ia;
        double ct = ia.getCollisionTime();
        if( ct < minCollisionTime) {
            minCollisionTime = ct;
            nextCollider = ia;
        }
    }
}

/**
 * Advances to next collision, applies collision dynamics to colliders and updates collision time/partners.
 * Also invokes collisionAction method of all registered integrator meters.
 */
protected void advanceToCollision() {
    advanceAcrossTimeStep(nextCollider.getCollisionTime());
    if(nextCollider.isFieldCollision) {  //one-body collision occurred
        Atom a = nextCollider.atom;
        nextCollider.collisionPotentialField.bump(a);
        upList(a);
        downList(a);
        
        //reset collision partners of atoms that are now up from this atom but still list it as their
        //collision partner.  Assumes this atom was moved down list, but this won't always be the case
        //This bit could be made more efficient
        upPairIterator.reset(a);
        while(upPairIterator.hasNext()) {
            AtomPair pair = upPairIterator.next();
            if(((Agent)pair.atom2().ia).collisionPartner == a) {  //upList atom could have atom as collision partner if atom was just moved down list
                upList(pair.atom2());
            }
        }
        //to keep collision lists perfect, should do an upList on atoms that had this
        //atom on its neighbor list, but no longer do because it has moved away
        
    }//end of if
    else {
        Atom partner = nextCollider.collisionPartner;
//        Atom partnerNextAtom = partner.nextMoleculeFirstAtom();  //put this back in for multiatomic speciesSwitch; also need to do more work with loop below
//        Atom partnerNextAtom = partner.coordinate.nextNeighbor().atom();
        Atom partnerNextAtom = null;  //remove this -- temporary
        atomPair.reset(nextCollider.atom,partner);
        Potential.Hard potential = nextCollider.collisionPotential;
        potential.bump(atomPair);
        for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
            cll.listener.collisionAction(atomPair, potential);
        }
        Atom a1N = atomPair.atom1();  //bump might have changed nextCollider or partner to new atoms in atomPair
        Atom a2P = atomPair.atom2();
//        nextCollider.getCollisionPotential().bump(nextCollider.atom,partner);
                
//        boolean upListedN = false;
//        boolean upListedP = false;
//        for(Atom a=firstPhase.firstAtom(); a!=partnerNextAtom; a=a.nextAtom()) {  //note that nextCollider's or partner's position in linked-list may have been moved by the bump method
//        upAtomIterator.reset(firstPhase.firstAtom());

//   Do upList for any atoms that were scheduled to collide with atoms colliding now
        upAtomIterator.reset();  //first atom in first cell
        while(upAtomIterator.hasNext()) {
            Atom a = upAtomIterator.next();
            if(a == a1N || a == a2P) continue;
            if(a == partnerNextAtom) break;
            Atom aPartner = ((Agent)a.ia).collisionPartner;
            if(aPartner==nextCollider.atom || aPartner==partner) {
                upList(a);
//                if(a == nextCollider.atom) {upListedN = true;}
//                else if(a == partner) {upListedP = true;}
            }
        }
//        if(!upListedN) {upList(atomPair.atom1());}  //nextCollider.atom
//        if(!upListedP) {upList(atomPair.atom2());}  //partner
        upList(a1N);
        upList(a2P);
        downList(a1N);
        downList(a2P);
    }//end of else

    findNextCollider();
}



/**
 * Advances all atom coordinates by tStep, without any intervening collisions.
 * Uses free-flight kinematics.
 */
protected void advanceAcrossTimeStep(double tStep) {
            
    for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
        ((Agent)a.ia).decrementCollisionTime(tStep);
        if(a.isStationary()) {continue;}  //skip if atom is stationary
        a.coordinate.freeFlight(tStep);
   //     a.translateBy(tStep*a.rm(),a.momentum());
    }
}

/**
 * Looks for collision of this atom with all atoms upList of the given atom, 
 * and sets this atom's agent to the next collision time and partner found there.
 */
protected void upList(Atom atom) {  
            
    double minCollisionTime = Double.MAX_VALUE;
    Agent aia = (Agent)atom.ia;
    aia.resetCollision();
            
    upPairIterator.reset(atom);
    while(upPairIterator.hasNext()) {
        AtomPair pair = upPairIterator.next();
        Potential.Hard potential = (Potential.Hard)parentSimulation().getPotential(pair);
        double time = potential.collisionTime(pair);
        if(time < minCollisionTime) {
            minCollisionTime = time;
            aia.setCollision(time,pair.atom2(),potential);  //atom2 is inner loop
        }
    }
      
    //Add in collisions with fields acting in the phase
    for(PotentialField f=firstPhase.firstField(); f!=null; f=f.nextField()) {
        if(f instanceof PotentialField.Hard) {
            PotentialField.Hard field = (PotentialField.Hard)f;
            double time = field.collisionTime(atom);
            if(time < minCollisionTime) {
                minCollisionTime = time;
                aia.setCollision(time,field);
            }
        }
    }
}

/**
 * Looks for collision of this atom with all atoms downList of the given atom, 
 * and updates the partner atom's agent if collision time with this is less than
 * its current collision time.
 */
protected void downList(Atom atom) {
            
    downPairIterator.reset(atom);
    while(downPairIterator.hasNext()) {
        AtomPair pair = downPairIterator.next();
        Agent aia = (Agent)pair.atom2().ia;  //atom2 is inner loop
        Potential.Hard potential = (Potential.Hard)parentSimulation().getPotential(pair);
        double time = potential.collisionTime(pair);
        if(time < aia.getCollisionTime()) {
            aia.setCollision(time,atom,potential);
        }
    }
}
          
/**
 * The simulation time elapsed since the start of the integration.
 * Cannot be reset to zero.  
 * Modified from superclass method to return correct value if in the middle
 * of processing collisions between time steps.
 */
public double elapsedTime() {return super.elapsedTime() + timeIncrement;}

/**
 * Sets up the integrator for action.
 * Do an upList call for each atom and find the next collider
 */
 protected void doReset() {
    if(isothermal) scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
    Atom.Iterator iterator = firstPhase.iteratorFactory().makeAtomIteratorUp();
    iterator.reset();
    while(iterator.hasNext()) {upList(iterator.next());}
    findNextCollider();
}
          
/**
 * Crude method to enforce constant-temperature constraint
 * Scales momenta of all atoms by a constant factor so that 
 * phase adheres to setpoint temperature.  Updates collision times appropriately.
 */
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
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
    }
    
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
    }

 /**
  * Produces the Agent defined by this integrator.
  * One instance of an Agent is placed in each atom controlled by this integrator.
  */
    public Integrator.Agent makeAgent(Atom a) {
        return new Agent(a);
    }
 
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
        public final double getCollisionTime() {return collisionTime;}
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
        public Color colliderColor = Color.red;
        /**
         * Color applied to the upList atom of the colliding pair
         */
        public Color partnerColor = Color.blue;
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
    
}

