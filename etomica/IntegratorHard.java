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
 * @see PotentialHard
 * @author David Kofke
 *
 */
 
public class IntegratorHard extends IntegratorHardAbstract implements EtomicaElement {

    public String getVersion() {return "IntegratorHard:01.06.14/"+IntegratorHardAbstract.VERSION;}
    
    Potential.Hard potential;

    public static final CollisionHandler collisionHandlerUp = new CollisionHandler() {
        double minCollisionTime;
        IntegratorHardAbstract.Agent aia;
        public void setAtom(Atom a) {
            minCollisionTime = Double.MAX_VALUE;
            aia = (IntegratorHardAbstract.Agent)a.ia;
            aia.resetCollision();
        }
        public void addCollision(Atom atom2, double collisionTime) {
            if(collisionTime < minCollisionTime) {
                minCollisionTime = collisionTime;
                aia.setCollision(collisionTime,atom2,potential);
            }
        }
    }; //end of collisionHandlerUp

    public static final CollisionHandler collisionHandlerDown = new CollisionHandler() {
        Atom atom1;
        public void setAtom(Atom a) {atom1 = a;}
        public void addCollision(Atom atom2, double collisionTime) {
            IntegratorHardAbstract.Agent aia = (IntegratorHardAbstract.Agent)atom2.ia;
            if(collisionTime < aia.collisionTime()) {
                aia.setCollision(collisionTime,atom1,potential);
            }
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
    
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return;
        potential = p.potential;
    }

	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
        upAtomIterator = factory.makeAtomIteratorUp();
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
                nextCollider = ia;
            }
        }
    }

    /**
    * Advances to next collision, applies collision dynamics to colliders and updates collision time/partners.
    * Also invokes collisionAction method of all registered integrator meters.
    */
    protected void processCollision() {
        if(nextCollider.isFieldCollision) {  //one-body collision occurred
            Atom a = nextCollider.atom;
            nextCollider.collisionPotentialField.bump(a);
            potential.findCollisions(a,   UP, collisionHandlerUp);
            potential.findCollisions(a, DOWN, collisionHandlerDown);
            
            //reset collision partners of atoms that are now up from this atom but still list it as their
            //collision partner.  Assumes this atom was moved down list, but this won't always be the case
            //This bit could be made more efficient
            
            //if(a movedInList) {  add a means for bump method to declare it moved atom in the list
            upAtomIterator.reset(a);
            while(upAtomIterator.hasNext()) {
                Atom atom = upAtomIterator.next();
                if(((Agent)atom.ia).collisionPartner == a) {  //upList atom could have atom as collision partner if atom was just moved down list
                    potential.findCollisions(atom, UP, collisionHandlerUp);
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
                    potential.findCollisions(a, UP, collisionHandlerUp);
    //                if(a == nextCollider.atom) {upListedN = true;}
    //                else if(a == partner) {upListedP = true;}
                }
            }
    //        if(!upListedN) {upList(atomPair.atom1());}  //nextCollider.atom
    //        if(!upListedP) {upList(atomPair.atom2());}  //partner

            potential.findCollisions(a1N,   UP, collisionHandlerUp);
            potential.findCollisions(a2P,   UP, collisionHandlerDown);
            potential.findCollisions(a1N, DOWN, collisionHandlerUp);
            potential.findCollisions(a2P, DOWN, collisionHandlerDown);

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
    * Sets up the integrator for action.
    * Do an upList call for each atom and find the next collider
    */
    protected void doReset() {
        if(isothermal) scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
        potential.findCollisions(collisionHandlerUp);
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

