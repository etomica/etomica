package etomica;

/**
 * General-purpose integrator for hard potentials.
 * Integrates equations of motion through time by advancing atoms from one collision to the
 * next.  Determination of time of collision and implementation of collision
 * dynamics is handled by the potential between the atoms, which therefore must
 * implement PotentialHard.  Each atom keeps in its Agent (obtained from this integrator) the shortest
 * time to collision with any of the atoms uplist of it (as defined by the iterator from 
 * the phase).
 *
 * @author David Kofke
 *
 */
public class IntegratorHard extends IntegratorHardAbstract implements EtomicaElement {

    public String getVersion() {return "IntegratorHard:01.06.27/"+IntegratorHardAbstract.VERSION;}
    
    private AtomIterator upAtomIterator;
    private static final IteratorDirective upList = new IteratorDirective(IteratorDirective.UP);
    private static final IteratorDirective downList = new IteratorDirective(IteratorDirective.DOWN);

    //collision handler is passed to the potential and is notified of each collision
    //the potential detects.  The collision handler contains the necessary logic to
    //process this information so that the collision lists are kept up to date.
    //The potential should call the handler's setPotential method, with itself as the
    //argument, before beginning to detect collisions. 
    
    //the up-handler has the logic of the Allen & Tildesley upList subroutine
    //sets collision time of given atom to minimum value for collisions with all atoms uplist of it
    public final IntegratorHardAbstract.CollisionHandler 
                    collisionHandlerUp = new IntegratorHardAbstract.CollisionHandler() {
        double minCollisionTime;
        IntegratorHardAbstract.Agent aia;
        Atom atom1;
        public IntegratorHardAbstract.CollisionHandler setAtom(Atom a) {
            atom1 = a;
            aia = (IntegratorHardAbstract.Agent)a.ia;
            minCollisionTime = aia.collisionTime();
            return this;
        }
        public void addCollision(AtomPair pair, double collisionTime) {
            if(pair.atom1() != atom1) setAtom(pair.atom1()); //need this if doing minimum collision time calculation for more than one atom
            if(collisionTime < minCollisionTime) {
                minCollisionTime = collisionTime;
                aia.setCollision(collisionTime, pair.atom2(), this.potential);
            }
        }
        public void addCollision(Atom atom, double collisionTime) {
            setAtom(atom);
            if(collisionTime < minCollisionTime) {
                minCollisionTime = collisionTime;
                aia.setCollision(collisionTime, null, this.potential);
            }
        }
    }; //end of collisionHandlerUp

    //the down-handler has the logic of the Allen & Tildesley downList subroutine
    //sets collision times of atoms downlist of given atom to minimum of their current
    //value and their value with given atom
    public final IntegratorHardAbstract.CollisionHandler 
                    collisionHandlerDown = new IntegratorHardAbstract.CollisionHandler() {
        public IntegratorHardAbstract.CollisionHandler setAtom(Atom a) {
            return this;
        }
        public void addCollision(AtomPair pair, double collisionTime) {
            IntegratorHardAbstract.Agent aia = (IntegratorHardAbstract.Agent)pair.atom2().ia;
            if(collisionTime < aia.collisionTime()) {
                aia.setCollision(collisionTime, pair.atom1(), this.potential);
            }
        }
        public void addCollision(Atom atom, double collisionTime) {//this shouldn't be called
            System.out.println("Unexpected entry into method addCollision in IntegratorHard's collisionHandlerDown");
        }
    }; //end of collisionHandlerDown

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
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactoryObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
	    super.makeIterators(factory);
        upAtomIterator = factory.makeAtomIterator();
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
            double ct = ia.collisionTime();
            if( ct < minCollisionTime) {
                minCollisionTime = ct;
                colliderAgent = ia;
            }
        }
    }

    /**
    * Updates collision times/partners for collider and partner, and 
    * for atoms that were to collide with one of them.
    */
    protected void updateCollisions() {
        
        Atom collider = colliderAgent.atom();
        Atom partner = colliderAgent.collisionPartner();
            
    //   Do upList for any atoms that were scheduled to collide with atoms colliding now
    //   Assumes collider and partner haven't moved in list
        upAtomIterator.reset();  //first atom in first cell
        while(upAtomIterator.hasNext()) {
            Atom a = upAtomIterator.next();
            if(a == collider) { //finished with atoms before collider...
                if(partner == null) break;  //and there is no partner, so we're done, or...
                else continue;              //...else just skip this atom and continue with loop
            }
            if(a == partner) break; //finished with atoms before partner; we're done
            IntegratorHardAbstract.Agent aAgent = (IntegratorHardAbstract.Agent)a.ia;
            Atom aPartner = aAgent.collisionPartner();
            if(aPartner == collider || aPartner == partner) {
                aAgent.resetCollision();
                phasePotential.findCollisions(upList.set(a), collisionHandlerUp.setAtom(a));
            }
        }//end while
            //reset collision partners of atoms that are now up from this atom but still list it as their
            //collision partner.  Assumes this atom was moved down list, but this won't always be the case
            //This bit could be made more efficient
            
            //if(a movedInList) {  add a means for bump method to declare it moved atom in the list
       /*     upAtomIterator.reset(a);
            while(upAtomIterator.hasNext()) {
                Atom atom = upAtomIterator.next();
                if(((Agent)atom.ia).collisionPartner == a) {  //upList atom could have atom as collision partner if atom was just moved down list
                    phasePotential.findCollisions(atom, UP, collisionHandlerUp);
                }
            }*/
            //to keep collision lists perfect, should do an upList on atoms that had this
            //atom on its neighbor list, but no longer do because it has moved away


        colliderAgent.resetCollision();
        phasePotential.findCollisions(upList.set(collider), collisionHandlerUp.setAtom(collider));
        phasePotential.findCollisions(downList.set(collider), collisionHandlerDown);
        if(partner != null) {
            ((IntegratorHardAbstract.Agent)partner.ia).resetCollision();
            phasePotential.findCollisions(upList.set(partner), collisionHandlerUp.setAtom(partner));
            phasePotential.findCollisions(downList.set(partner), collisionHandlerDown);
        }

    }//end of processCollision

    /**
    * Advances all atom coordinates by tStep, without any intervening collisions.
    * Uses free-flight kinematics.
    */
    protected void advanceAcrossTimeStep(double tStep) {
        
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
     //   for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            a.coordinate.freeFlight(tStep);
    //     a.translateBy(tStep*a.rm(),a.momentum());
        }
    }

    /**
    * Sets up the integrator for action.
    * Do an upList call for each atom and find the next collider
    */
    protected void doReset() {
        if(isothermal) scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            ((IntegratorHardAbstract.Agent)atomIterator.next().ia).resetCollision();
        }
        phasePotential.findCollisions(upList.set(), collisionHandlerUp);
        findNextCollider();
    }
                
    /**
    * Produces the Agent defined by this integrator.
    * One instance of an Agent is placed in each atom controlled by this integrator.
    */
        public Integrator.Agent makeAgent(Atom a) {
            return new IntegratorHardAbstract.Agent(a);
        }
             
}//end of IntegratorHard

