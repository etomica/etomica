package etomica.integrator;

import etomica.Atom;
import etomica.AtomsetIterator;
import etomica.Debug;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.atom.AtomArrayList;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialHard;
import etomica.space.CoordinatePair;
import etomica.space.ICoordinateKinetic;
import etomica.space.Vector;
import etomica.utility.TreeLinker;
import etomica.utility.TreeList;

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
 
public class IntegratorHard extends IntegratorMD {

    //handle to the integrator agent holding information about the next collision
    protected IntegratorHard.Agent colliderAgent;
    //first of a linked list of objects (typically meters) that are called each time a collision is processed
    protected CollisionListenerLinker collisionListenerHead = null;
    private Atom[] atoms;
    Vector c3;
    CoordinatePair cPairDebug;

    protected final IteratorDirective upList = new IteratorDirective(IteratorDirective.UP);
    protected final IteratorDirective downList = new IteratorDirective(IteratorDirective.DOWN);
    protected final ReverseCollisionHandler reverseCollisionHandler = new ReverseCollisionHandler();

    protected final CollisionHandlerUp collisionHandlerUp = new CollisionHandlerUp();
    protected final CollisionHandlerDown collisionHandlerDown = new CollisionHandlerDown();
    protected final AtomArrayList listToUpdate = new AtomArrayList();
    
    protected PotentialHard nullPotential;
    protected double collisionTimeStep;
    protected int collisionCount;

    protected TreeList eventList;

    private Atom[] targetAtom = new Atom[1];
    
    public IntegratorHard(PotentialMaster potentialMaster) {
        super(potentialMaster);
        nullPotential = null;
        atoms = new Atom[2];
        eventList = new TreeList();
    }
    
    public boolean addPhase(Phase phase) {
        if(!super.addPhase(phase)) return false;
        atomIterator.setPhase(phase);
        return true;
    }
    
    public IntegratorHard.Agent colliderAgent() {
        return colliderAgent;
    }
    
    /** 
     * Steps all atoms across time interval timeStep, handling all intervening collisions.
     */
    public void doStep() {
        findNextCollider();
        collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.MAX_VALUE;
        double oldTime = 0;
        while(collisionTimeStep < timeStep) {//advance to collision if occurs before remaining interval
            atoms[0] = colliderAgent.atom();
            atoms[1] = colliderAgent.collisionPartner();
            if (collisionTimeStep - oldTime < -1.e-10/Math.sqrt(temperature)) {
                System.out.println("previous collision occured before current one");
                System.out.println("previous time: "+oldTime+" current time: "+collisionTimeStep);
                System.out.println("collision between "+atoms[0]+" and "+atoms[1]+" potential "+colliderAgent.collisionPotential.getClass());
                if (atoms[1] != null) {
                    cPairDebug = Simulation.getDefault().space.makeCoordinatePair();
                    cPairDebug.setNearestImageTransformer(firstPhase.boundary());
//                    cPairDebug.trueReset(atoms[0].coord,atoms[1].coord,oldTime);
                    System.out.println("distance at last collision time was "+Math.sqrt(cPairDebug.r2()));
//                    cPairDebug.trueReset(collisionTimeStep);
                    System.out.println("distance now "+Math.sqrt(cPairDebug.r2()));
                }
                throw new RuntimeException("this simulation is not a time machine");
            }
            if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.thisPhase(firstPhase)) || Debug.anyAtom(atoms))) {
                System.out.println("collision between atoms "+atoms[0]+" and "+atoms[1]+" at "+collisionTimeStep);
            }
            if (Debug.ON && Debug.DEBUG_NOW && Debug.ATOM1 != null 
                  && Debug.ATOM2 != null && Debug.thisPhase(firstPhase)) {
                cPairDebug = Simulation.getDefault().space.makeCoordinatePair();
                cPairDebug.setNearestImageTransformer(firstPhase.boundary());
//                cPairDebug.trueReset(Debug.ATOM1.coord,Debug.ATOM2.coord,collisionTimeStep);
                double r2 = cPairDebug.r2();
                if (Debug.LEVEL > 1 || Math.sqrt(r2) < Default.ATOM_SIZE-1.e-11) {
                    System.out.println("distance between "+Debug.ATOM1+" and "+Debug.ATOM2+" is "+Math.sqrt(r2));
                    if (Debug.LEVEL > 2 || Math.sqrt(r2) < Default.ATOM_SIZE-1.e-11) {
//                        System.out.println(Debug.ATOM1+" coordinates "+Debug.ATOM1.coord.truePosition(collisionTimeStep));
//                        System.out.println(Debug.ATOM2+" coordinates "+Debug.ATOM2.coord.truePosition(collisionTimeStep));
                    }
                }
                Debug.checkAtoms();
                if (Debug.LEVEL > 1) {
                    PotentialHard p = ((Agent)Debug.ATOM1.ia).collisionPotential;
                    System.out.println(Debug.ATOM1+" collision time "+((Agent)Debug.ATOM1.ia).collisionTime+" with "+((Agent)Debug.ATOM1.ia).collisionPartner
                            +" with potential "+(p!=null ? p.getClass() : null));
                }
            }

            colliderAgent.collisionPotential.bump(atoms,collisionTimeStep);
            double dE = colliderAgent.collisionPotential.energyChange();
            currentPotentialEnergy[0] += dE;
            currentKineticEnergy[0] -= dE;
            
            for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
                cll.listener.collisionAction(colliderAgent);
            }
            collisionCount++;
            if (atoms[1] == null) {
                updateAtom(atoms[0]);
            }
            else {
                updateAtoms(atoms);
            }
            findNextCollider(); //this sets colliderAgent for the next collision
            
            oldTime = collisionTimeStep;
            collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.MAX_VALUE;
        } 
        advanceAcrossTimeStep(timeStep);
        collisionTimeStep = 0.0;
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.thisPhase(firstPhase)) {
            eventList.check();
            meterPE.setPhase(phase);
            double[] PE = meterPE.getData();
            for (int i=0; i<phase.length; i++) {
                if (PE[i] != currentPotentialEnergy[i]) {
                    throw new RuntimeException("potential energy of "+phase[i]+" is wrong. it's actually "+PE[i]+" but I thought it was "+currentPotentialEnergy[i]);
                }
            }
            meterKE.setPhase(phase);
            double[] KE = meterKE.getData();
            for (int i=0; i<phase.length; i++) {
                if (Math.abs(KE[i] - currentKineticEnergy[i]) > 1.e-9) {
                    throw new RuntimeException("kinetic energy of "+phase[i]+" is wrong. it's actually "+KE[i]+" but I thought it was "+currentKineticEnergy[i]);
                }
            }
        }

        if(isothermal) doThermostat();
    }//end of doStep
    
    public int getCollisionCount() {
        return collisionCount;
    }
    
   /**
	* Loops through all atoms to identify the one with the smallest value of collisionTime
	* Collision time is obtained from the value stored in the Integrator.Agent from each atom.
	*/
	protected void findNextCollider() {
        colliderAgent = (Agent)eventList.firstElement();
	}

    /**
     * Updates collision times/partners for collider and partner, and 
     * for atoms that were to collide with one of them.  This method 
     * should not be called for 1-body collisions (call updateAtom).
     */
    protected void updateAtoms(Atom[] colliders) {

        listToUpdate.clear();

        targetAtom[0] = colliders[0];
        downList.setTargetAtoms(targetAtom);
        potential.calculate(firstPhase, downList, reverseCollisionHandler);
        targetAtom[0] = colliders[1];
        downList.setTargetAtoms(targetAtom);
        potential.calculate(firstPhase, downList, reverseCollisionHandler);

        // this will also update colliders[0] since it thought it would collide 
        // with colliders[1] (and did)
        processReverseList();

        targetAtom[0] = colliders[0];
        downList.setTargetAtoms(targetAtom);
        potential.calculate(firstPhase, downList, collisionHandlerDown);

        Agent agent = (Agent)colliders[1].ia;
        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollision();
        targetAtom[0] = colliders[1];
        upList.setTargetAtoms(targetAtom);
        potential.calculate(firstPhase, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        downList.setTargetAtoms(targetAtom);
        collisionHandlerUp.setAtom(colliders[1]);
        potential.calculate(firstPhase, downList, collisionHandlerDown);

    }
    

    /**
     * updates collision time for a single atom (and any atom that might a 
     * collision with the atom)
     */
    protected void updateAtom(Atom a) {

        listToUpdate.clear();

        targetAtom[0] = a;
        downList.setTargetAtoms(targetAtom);
        potential.calculate(firstPhase, downList, reverseCollisionHandler);
        processReverseList();

        Agent agent = (Agent)a.ia;
        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollision();
        targetAtom[0] = a;
        upList.setTargetAtoms(targetAtom);
        downList.setTargetAtoms(targetAtom);
        collisionHandlerUp.setAtom(a);
        potential.calculate(firstPhase, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        potential.calculate(firstPhase, downList, collisionHandlerDown);
    }

    /**
     * Recalculate collision times (up) for any atoms atoms in the "reverse"
     * list of atoms that thought they would collide with an atom that's
     * changed.
     */
    private void processReverseList() {
        int size = listToUpdate.size();
        for (int i=0; i<size; i++) {
            Atom reverseAtom = listToUpdate.get(i);
            Agent agent = (Agent)reverseAtom.ia;
            if (agent.collisionPotential != null) {
                agent.eventLinker.remove();
            }
            agent.resetCollision();
            targetAtom[0] = reverseAtom;
            upList.setTargetAtoms(targetAtom);
            collisionHandlerUp.setAtom(reverseAtom);
            potential.calculate(firstPhase, upList, collisionHandlerUp);
            if (agent.collisionPotential != null) {
                eventList.add(agent.eventLinker);
            }
        }
    }
    

	/**
     * Advances all atom coordinates by tStep, without any intervening collisions.
     * Uses free-flight kinematics.
     */
	protected void advanceAcrossTimeStep(double tStep) {
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			Atom a = atomIterator.nextAtom();
			((Agent)a.ia).decrementCollisionTime(tStep);
			a.coord.position().PEa1Tv1(tStep,((ICoordinateKinetic)a.coord).velocity());
		}
	}

    public void reset() {
        super.reset();
        neighborsUpdated();
    }

    /**
     * Do an upList call for each atom and reconstruct the event list.
     */
    public void neighborsUpdated() {
        super.neighborsUpdated();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            ((Agent)atomIterator.nextAtom().ia).resetCollision();
        }
        targetAtom[0] = null;
        upList.setTargetAtoms(targetAtom);
        collisionHandlerUp.reset();
        potential.calculate(firstPhase, upList, collisionHandlerUp); //assumes only one phase
        eventList.reset();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Agent agent = (Agent)atomIterator.nextAtom().ia;
            if (agent.collisionPotential != null) {
                eventList.add(agent.eventLinker);
            }
        }
    }

    /**
     * Updates collision times appropriately after scaling momenta.
     */
    protected double scaleMomenta(Phase aPhase) {
        double s = super.scaleMomenta(aPhase);
        double rs = 1.0/s;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            ((Agent)a.ia).collisionTime *= rs;
            ((Agent)a.ia).eventLinker.sortKey = ((Agent)a.ia).collisionTime; 
        }
        // don't need to update eventTree because event order didn't change
        return s;
    }

    /**
     * Updates collision times appropriately after randomizing momenta
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomenta(Phase aPhase) {
        super.randomizeMomenta(aPhase);
        // super.randomizeMomenta recalculates the kinetic energy and doesn't
        // change the potential energy, so just act like the neighbor lists were
        // updated
        neighborsUpdated();
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
     * @return Returns the nullPotential.
     */
    public PotentialHard getNullPotential() {
        return nullPotential;
    }
    /**
     * sets the "null" potential, a 1-body potential used
     * by IntegratorHard to prevent Atoms from wrapping around
     * periodic boundaries when neighbor listing is not used.
     * @param nullPotential The nullPotential to set.
     */
    public void setNullPotential(PotentialHard nullPotential) {
        this.nullPotential = nullPotential;
    }
    
    //collision handler is passed to the potential and is notified of each collision
    //the potential detects.  The collision handler contains the necessary logic to
    //process this information so that the collision lists are kept up to date.
    //The potential should call the handler's setPotential method, with itself as the
    //argument, before beginning to detect collisions. 

    //the up-handler has the logic of the Allen & Tildesley upList subroutine
    //sets collision time of given atom to minimum value for collisions with all atoms uplist of it
    private final class CollisionHandlerUp extends PotentialCalculation {
        double minCollisionTime;
        IntegratorHard.Agent aia;
        Atom atom1;

        /**
         * resets the "atom" held by this class, ensuring the method will
         * work even if (coincidentally) the atom is the same as the same 
         * one as the last time the method was called.
         */
        public void reset() {
            atom1 = null;
        }

        /**
         * sets the atom whose collision time is to be calculated
         */
        public void setAtom(Atom a) {
            atom1 = a;
            aia = (IntegratorHard.Agent)a.ia;
            minCollisionTime = aia.collisionTime();
        }//end of setAtom

        //atom pair
        public void doCalculation(AtomsetIterator iterator, Potential potential) {
            iterator.reset();
            PotentialHard pHard = (PotentialHard)potential; 
            while (iterator.hasNext()) {
                Atom[] atoms = iterator.next();
                if(atoms[0] != atom1) setAtom(atoms[0]); //need this if doing minimum collision time calculation for more than one atom
                double collisionTime = pHard.collisionTime(atoms,collisionTimeStep);
                if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atoms)) || Debug.allAtoms(atoms))) {
                    System.out.println("collision up time "+collisionTime+" for atom "+atoms[0]+(atoms.length > 1 ? " with "+atoms[1] : ""));
                }
                if(collisionTime < minCollisionTime) {
                    if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atoms))) {
                        System.out.println("setting up time "+collisionTime+" for atom "+atoms[0]+(atoms.length > 1 ? " with "+atoms[1] : ""));
                    }
                    minCollisionTime = collisionTime;
                    aia.setCollision(collisionTime, atoms.length == 1 ? null : atoms[1], pHard);
                }//end if
            }
        }//end of calculate(AtomPair...

    } //end of collisionHandlerUp

	//the down-handler has the logic of the Allen & Tildesley downList subroutine
	//sets collision times of atoms downlist of given atom to minimum of their current
	//value and their value with given atom
	private final class CollisionHandlerDown extends PotentialCalculation {
		public void doCalculation(AtomsetIterator iterator, Potential potential) {
			if (potential.nBody() != 2) return;
			iterator.reset();
            PotentialHard pHard = (PotentialHard)potential;
			while (iterator.hasNext()) {
				Atom[] atomPair = iterator.next();
				double collisionTime = pHard.collisionTime(atomPair,collisionTimeStep);
				if(collisionTime < Double.MAX_VALUE) {
					Agent aia = (Agent)atomPair[1].ia;
					if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atomPair)))) {
						System.out.println("collision down time "+collisionTime+" for atom "+atomPair[1]+" with "+atomPair[0]);
					}
					if(collisionTime < aia.collisionTime()) {
						if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atomPair))) {
							System.out.println("setting down time "+collisionTime+" for atom "+atomPair[1]+" with "+atomPair[0]);
						}
                        if (aia.collisionPotential != null) {
                            aia.eventLinker.remove();
                        }
                        aia.setCollision(collisionTime, atomPair[0], pHard);
                        eventList.add(aia.eventLinker);
					}//end if
				}//end if
			}
		}//end of actionPerformed
	} //end of collisionHandlerDown


    /**
     * A PotentialCalculation to find atoms that thought they would collide
     * with an atom.  The iterator should return an atom and its "down"
     * neighbors.  
     */
    private final class ReverseCollisionHandler extends PotentialCalculation {
        
        public void doCalculation(AtomsetIterator iterator, Potential p) {
            if (iterator.nBody() != 2) return;
            iterator.reset();
            // look for pairs in which pair[0] is the collision partner of pair[1]
            while (iterator.hasNext()) {
                Atom[] pair = iterator.next();
                Agent aAgent = (Agent)pair[1].ia;
                Atom aPartner = aAgent.collisionPartner();
                if (Debug.ON && Debug.DEBUG_NOW && ((Debug.allAtoms(pair) && Debug.LEVEL > 1) || (Debug.anyAtom(pair) && Debug.LEVEL > 2))) {
                    System.out.println(pair[1]+" thought it would collide with "+aPartner);
                }
                if(aPartner == pair[0]) {
                    if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(pair) || Debug.LEVEL > 1)) {
                        System.out.println("Will update "+pair[1]+" because it wanted to collide with "+aPartner);
                    }
//                    Evil.reverseAtomLinker++;
                    listToUpdate.add(pair[1]);
                }
            }
        }
    }

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("Collision-based molecular dynamics simulation of hard potentials");
		return info;
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
    }//end of CollisionListenerLinker

	/**
	* Produces the Agent defined by this integrator.
	* One instance of an Agent is placed in each atom controlled by this integrator.
	*/
		public Object makeAgent(Atom a) {
			return new IntegratorHard.Agent(a);
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
    public class Agent {  //need public so to use with instanceof
        public Atom atom, collisionPartner;
        public double collisionTime = Double.MAX_VALUE; //time to next collision
        public PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        public TreeLinker eventLinker;
        
        public Agent(Atom a) {
            atom = a;
            eventLinker = new TreeLinker(this);
        }
        
        public String toString() {
            return "Collider: "+atom.toString()+"; Partner: "+collisionPartner.toString()+"; Potential: "+collisionPotential.toString();
        }
        public final Atom atom() {return atom;}
        public final Atom collisionPartner() {return collisionPartner;}

        /**
         * resets time, potential and partner.  caller should remove 
         * eventLinker from the tree if needed before calling this method.
         */
        public void resetCollision() {
            // events with a (non-null) potential must be in the tree
            // events in the tree must have a potential
            collisionPotential = nullPotential;
            if (nullPotential != null) {
                targetAtom[0] = atom;
                collisionTime = nullPotential.collisionTime(targetAtom,collisionTimeStep);
            }
            else {
                collisionTime = Double.MAX_VALUE;
            }
            collisionPartner = null;
            eventLinker.sortKey = collisionTime;
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
            eventLinker.sortKey = time;
        }


    /**
     * Decreases the recorded time to collision of this atom
     * This action is performed when the atom is advanced without a collision
     */
        public final void decrementCollisionTime(double interval) {
            collisionTime -= interval;
            eventLinker.sortKey = collisionTime;
//            time0 += interval;
        }
    /**
     * Accessor method for the time to next collision of this atom
     */
        public final double collisionTime() {return collisionTime;}
    }//end of Agent
    
    public interface CollisionListener {
        public void collisionAction(Agent colliderAgent);
    }

    
	/**
	 * Demonstrates how this class is implemented.
	 */
 /*   public static void main(String[] args) {
		Simulation.instance = new etomica.graphics.SimulationGraphic(new Space2D());
		IntegratorHard integratorHard1 = new IntegratorHard();
	 //   integratorHard1.setTimeStep(0.02);
		SpeciesSpheresMono speciesSpheres1 = new SpeciesSpheresMono(10);
		SpeciesSpheresMono speciesSpheres2 = new SpeciesSpheresMono(1);
		etomica.graphics.ColorSchemeByType.setColor(speciesSpheres2, java.awt.Color.red);
		final Phase phase = new Phase();
		P2HardSphere potential12 = new P2HardSphere(4.0);
		P2HardSphere potential22 = new P2HardSphere(5.0);
		P2HardSphere potential11 = new P2HardSphere(3.0);
		speciesSpheres2.setDiameter(5.0);
		Controller controller1 = new Controller();
		etomica.graphics.DisplayPhase displayPhase1 = new etomica.graphics.DisplayPhase();
		etomica.graphics.DisplayTimer timer = new etomica.graphics.DisplayTimer(integratorHard1);
		timer.setUpdateInterval(10);
		MeterRDF meterRDF = new MeterRDF();
		etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot();
//		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go(); 
		
		potential12.setSpecies(speciesSpheres1, speciesSpheres2);
		potential22.setSpecies(speciesSpheres2, speciesSpheres2);
		potential11.setSpecies(speciesSpheres1, speciesSpheres1);
		/*
		potential.setIterator(new AtomPairIterator(Simulation.instance.space,
				new AtomIteratorSequential(speciesSpheres1.getAgent(phase),true),
				new AtomIteratorSequential(speciesSpheres2.getAgent(phase),true)));
                
		potential2.setIterator(new AtomPairIterator(Simulation.instance.space,
				new AtomIteratorSequential(speciesSpheres2.getAgent(phase),true),
				new AtomIteratorSequential(speciesSpheres2.getAgent(phase),true)));
                
		potential0.setIterator(new AtomPairIterator(Simulation.instance.space,
				new AtomIteratorSequential(speciesSpheres1.getAgent(phase),true),
				new AtomIteratorSequential(speciesSpheres1.getAgent(phase),true)));
				* /
	//    displayPhase1.setColorScheme(integratorHard1.new HighlightColliders());
	    
		displayPhase1.addDrawable(new etomica.graphics.Drawable() {
			public void draw(java.awt.Graphics g, int[] origin, double s) {
				double toPixels = etomica.units.BaseUnit.Length.Sim.TO_PIXELS*s;
				int i=0;
				for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {
					IntegratorHard.Agent agent = (IntegratorHard.Agent)a.ia;
					if(agent == null) return;
					String text = Float.toString((float)agent.collisionTime);
					Space.Vector r = a.coord.position();
					int xP = origin[0] + (int)(toPixels*(r.component(0)));
					int yP = origin[1] + (int)(toPixels*(r.component(1)));
					g.setColor(java.awt.Color.gray);
					g.drawString(text, xP, yP-20);
					g.setColor(java.awt.Color.red);
					g.drawString(Integer.toString(a.node.index()), xP-20, yP-20);
				}
			}
		});
		/*
		//writes collision data to console
		integratorHard1.addCollisionListener(new CollisionListener() {
			public void collisionAction(IntegratorHard.Agent agent) {
				int i0 = agent.atom.debugIndex;
				int i1 = (agent.collisionPartner!=null)?agent.collisionPartner.debugIndex:-1;
	  //          System.out.print(i0+" "+i1+" ");
				for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {
					IntegratorHard.Agent aAgent = (IntegratorHard.Agent)a.ia;
					if(aAgent == null) return;
					String time = Float.toString((float)aAgent.collisionTime);
					int idx = (aAgent.collisionPartner!=null)?aAgent.collisionPartner.debugIndex:-1;
					System.out.print("("+a.debugIndex+","+idx+","+time+") ");
				}
				System.out.println();
			}
		});* /
		etomica.graphics.SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
	}//end of main
*/
}//end of IntegratorHard

