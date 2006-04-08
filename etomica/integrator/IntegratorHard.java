package etomica.integrator;

import java.io.Serializable;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.Potential1;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialHard;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.TreeLinker;
import etomica.util.TreeList;

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
 
public class IntegratorHard extends IntegratorMD implements AgentSource {

    //handle to the integrator agent holding information about the next collision
    protected IntegratorHard.Agent colliderAgent;
    //first of a linked list of objects (typically meters) that are called each time a collision is processed
    protected CollisionListenerLinker collisionListenerHead = null;
    private final AtomPair pair;
    private double minDelta;
    private AtomPair debugPair;

    protected final IteratorDirective upList = new IteratorDirective(IteratorDirective.Direction.UP);
    protected final IteratorDirective downList = new IteratorDirective(IteratorDirective.Direction.DOWN);
    protected final AtomArrayList listToUpdate = new AtomArrayList();
    protected final TreeList eventList = new TreeList();

    protected final ReverseCollisionHandler reverseCollisionHandler;
    protected final CollisionHandlerUp collisionHandlerUp;
    protected final CollisionHandlerDown collisionHandlerDown;
    
    protected PotentialHard nullPotential;
    protected double collisionTimeStep;
    protected int collisionCount;
    
    protected Agent[] agents;
    protected AtomAgentManager agentManager;

    public IntegratorHard(Simulation sim) {
        this(sim.potentialMaster,sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorHard(PotentialMaster potentialMaster, double timeStep, 
            double temperature) {
        super(potentialMaster,timeStep,temperature);
        nullPotential = null;
        pair = new AtomPair();
        agents = new Agent[0];
        reverseCollisionHandler = new ReverseCollisionHandler(listToUpdate);
        reverseCollisionHandler.setAgents(agents);
        collisionHandlerUp = new CollisionHandlerUp();
        collisionHandlerUp.setAgents(agents);
        collisionHandlerDown = new CollisionHandlerDown(eventList);
        collisionHandlerDown.setAgents(agents);
    }
    
    public void setPhase(Phase newPhase) {
        if (phase != null) {
            // allow agentManager to de-register itself as a PhaseListener
            agentManager.setPhase(null);
        }
        super.setPhase(newPhase);
        agentManager = new AtomAgentManager(this,newPhase);
        agents = (Agent[])agentManager.getAgents();
        if (nullPotential != null) {
            ((Potential1)nullPotential).setPhase(newPhase);
        }
    }
    
    public void setup() {
        agents = (Agent[])agentManager.getAgents();
        super.setup();
    }
    
    /* (non-Javadoc)
     * @see etomica.Integrator#setTemperature(double)
     */
    public void setTemperature(double temperature) {
        super.setTemperature(temperature);
        minDelta = -5.e-8/Math.sqrt(temperature);
    }
    
    public IntegratorHard.Agent colliderAgent() {
        return colliderAgent;
    }
    
    /** 
     * Steps all atoms across time interval timeStep, handling all intervening collisions.
     */
    public void doStep() {
        findNextCollider();
        collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.POSITIVE_INFINITY;
        double oldTime = 0;
        while(collisionTimeStep < timeStep) {//advance to collision if occurs before remaining interval
            AtomSet atoms;
            if (colliderAgent.collisionPartner() != null) {
                atoms = pair;
                pair.atom0 = colliderAgent.atom();
                pair.atom1 = colliderAgent.collisionPartner();
            }
            else {
                atoms = colliderAgent.atom();
            }
            if (collisionTimeStep - oldTime < minDelta) {
                System.out.println("diff "+(collisionTimeStep - oldTime)+" minDelta "+minDelta);
                System.out.println("previous collision occured before current one");
                System.out.println("previous time: "+oldTime+" current time: "+collisionTimeStep);
                System.out.println("collision for "+atoms+" potential "+colliderAgent.collisionPotential.getClass());
                if (atoms instanceof AtomPair) {
                    Vector dr = phase.space().makeVector();
                    Vector dv = phase.space().makeVector();

                    AtomLeaf atom0 = (AtomLeaf)pair.atom0;
                    AtomLeaf atom1 = (AtomLeaf)pair.atom1;
                    ICoordinateKinetic coord0 = (ICoordinateKinetic)atom0.coord;
                    ICoordinateKinetic coord1 = (ICoordinateKinetic)atom1.coord;
                    dv.Ev1Mv2(coord1.velocity(), coord0.velocity());
                    
                    dr.Ev1Mv2(coord1.position(), coord0.position());
                    phase.getBoundary().nearestImage(dr);

                    dr.PEa1Tv1(oldTime,dv);
                    System.out.println("distance at last collision time was "+Math.sqrt(dr.squared()));
                    dr.PEa1Tv1(collisionTimeStep-oldTime,dv);
                    System.out.println("distance now "+Math.sqrt(dr.squared()));
                }
                throw new RuntimeException("this simulation is not a time machine");
            }
            if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.thisPhase(phase)) || Debug.anyAtom(atoms))) {
                System.out.println("collision between atoms "+atoms+" at "+collisionTimeStep+" with "+colliderAgent.collisionPotential.getClass());
            }
            if (Debug.ON && Debug.DEBUG_NOW && Debug.thisPhase(phase)) {
                debugPair = Debug.getAtoms(phase);
                if (debugPair.atom0 instanceof AtomLeaf && debugPair.atom1 instanceof AtomLeaf) {
                    Vector dr = phase.space().makeVector();
                    Vector dv = phase.space().makeVector();

                    AtomLeaf atom0 = (AtomLeaf)pair.atom0;
                    AtomLeaf atom1 = (AtomLeaf)pair.atom1;
                    ICoordinateKinetic coord0 = (ICoordinateKinetic)atom0.coord;
                    ICoordinateKinetic coord1 = (ICoordinateKinetic)atom1.coord;
                    dv.Ev1Mv2(coord1.velocity(), coord0.velocity());
                    
                    dr.Ev1Mv2(coord1.position(), coord0.position());
                    phase.getBoundary().nearestImage(dr);

                    dr.PEa1Tv1(collisionTimeStep,dv);
                    double r2 = dr.squared();
                    if (Debug.LEVEL > 1 || Math.sqrt(r2) < Debug.ATOM_SIZE-1.e-11) {
                        System.out.println("distance between "+debugPair+" is "+Math.sqrt(r2));
                        if (Debug.LEVEL > 2 || Math.sqrt(r2) < Debug.ATOM_SIZE-1.e-11) {
                            dr.Ea1Tv1(collisionTimeStep,((ICoordinateKinetic)((AtomLeaf)debugPair.atom0).coord).velocity());
                            System.out.println(debugPair.atom0+" coordinates "+((AtomLeaf)debugPair.atom0).coord.position().P(dr));
                            dr.Ea1Tv1(collisionTimeStep,((ICoordinateKinetic)((AtomLeaf)debugPair.atom1).coord).velocity());
                            System.out.println(debugPair.atom1+" coordinates "+((AtomLeaf)debugPair.atom1).coord.position().P(dr));
                        }
                    }
                    if (Debug.LEVEL > 1) {
                        PotentialHard p = agents[debugPair.atom0.getGlobalIndex()].collisionPotential;
                        System.out.println(debugPair.atom0+" collision time "+agents[debugPair.atom0.getGlobalIndex()].collisionTime()+" with "+agents[debugPair.atom0.getGlobalIndex()].collisionPartner
                                +" with potential "+(p!=null ? p.getClass() : null));
                    }
                }
                else if (Debug.LEVEL > 2 && debugPair.atom0 instanceof AtomLeaf) {
                    Vector dr = potential.getSpace().makeVector();
                    dr.Ea1Tv1(collisionTimeStep,((ICoordinateKinetic)((AtomLeaf)debugPair.atom0).coord).velocity());
                    System.out.println(debugPair.atom0+" coordinates "+((AtomLeaf)debugPair.atom0).coord.position().P(dr));
                }
            }

            colliderAgent.collisionPotential.bump(atoms,collisionTimeStep);
            double dE = colliderAgent.collisionPotential.energyChange();
            currentPotentialEnergy += dE;
            currentKineticEnergy -= dE;
            
            for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
                cll.listener.collisionAction(colliderAgent);
            }
            collisionCount++;
            if (atoms instanceof Atom) {
                updateAtom((Atom)atoms);
            }
            else {
                updateAtoms((AtomPair)atoms);
            }
            findNextCollider(); //this sets colliderAgent for the next collision
            
            oldTime = collisionTimeStep;
            collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.POSITIVE_INFINITY;
        } 
        colliderAgent = null;

        advanceAcrossTimeStep(timeStep);
        collisionTimeStep = 0.0;
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.thisPhase(phase)) {
            eventList.check();
            meterPE.setPhase(phase);
            double PE = meterPE.getDataAsScalar();
            if (Math.abs((PE - currentPotentialEnergy)/(PE+currentPotentialEnergy)) > 1.e-9
                    && Math.abs(PE - currentPotentialEnergy) > 1.e-9) {
                throw new RuntimeException("potential energy of "+phase+" is wrong. it's actually "+PE+" but I thought it was "+currentPotentialEnergy);
            }
            meterKE.setPhase(phase);
            double KE = meterKE.getDataAsScalar();
            if (Math.abs(KE - currentKineticEnergy) > 1.e-8) {
                throw new RuntimeException("kinetic energy of "+phase+" is wrong. it's actually "+KE+" but I thought it was "+currentKineticEnergy);
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
    protected void updateAtoms(AtomPair colliders) {

        listToUpdate.clear();

        downList.setTargetAtom(colliders.atom0);
        potential.calculate(phase, downList, reverseCollisionHandler);
        downList.setTargetAtom(colliders.atom1);
        potential.calculate(phase, downList, reverseCollisionHandler);

        // this will also update colliders[0] since it thought it would collide 
        // with colliders[1] (and did)
        processReverseList();

        downList.setTargetAtom(colliders.atom0);
        collisionHandlerDown.collisionTimeStep = this.collisionTimeStep;
        potential.calculate(phase, downList, collisionHandlerDown);

        Agent agent = agents[colliders.atom1.getGlobalIndex()];
        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollision();
        upList.setTargetAtom(colliders.atom1);
        collisionHandlerUp.setAtom(colliders.atom1);
        collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
        potential.calculate(phase, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        downList.setTargetAtom(colliders.atom1);
        collisionHandlerDown.collisionTimeStep = this.collisionTimeStep;
        potential.calculate(phase, downList, collisionHandlerDown);

    }
    

    /**
     * updates collision time for a single atom (and any atom that might a 
     * collision with the atom)
     */
    protected void updateAtom(Atom a) {

        listToUpdate.clear();

        downList.setTargetAtom(a);
        potential.calculate(phase, downList, reverseCollisionHandler);
        processReverseList();

        Agent agent = agents[a.getGlobalIndex()];
        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollision();
        upList.setTargetAtom(a);
        collisionHandlerUp.setAtom(a);
        collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
        potential.calculate(phase, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        collisionHandlerDown.collisionTimeStep = this.collisionTimeStep;
        potential.calculate(phase, downList, collisionHandlerDown);
    }

    /**
     * Recalculate collision times (up) for any atoms atoms in the "reverse"
     * list of atoms that thought they would collide with an atom that's
     * changed.
     */
    protected void processReverseList() {
        int size = listToUpdate.size();
        for (int i=0; i<size; i++) {
            Atom reverseAtom = listToUpdate.get(i);
            Agent agent = agents[reverseAtom.getGlobalIndex()];
            if (agent.collisionPotential != null) {
                agent.eventLinker.remove();
            }
            agent.resetCollision();
            upList.setTargetAtom(reverseAtom);
            collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
            collisionHandlerUp.setAtom(reverseAtom);
            potential.calculate(phase, upList, collisionHandlerUp);
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
			AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            agents[a.getGlobalIndex()].decrementCollisionTime(tStep);
			a.coord.position().PEa1Tv1(tStep,((ICoordinateKinetic)a.coord).velocity());
		}
	}

    public void reset() throws ConfigurationOverlapException {
        // reset might be called because atoms were added or removed
        // calling getAgents ensures we have an up-to-date array.
        agents = (Agent[])agentManager.getAgents();
        collisionHandlerUp.setAgents(agents);
        collisionHandlerDown.setAgents(agents);
        reverseCollisionHandler.setAgents(agents);

        ConfigurationOverlapException overlapException = null;
        try {
            super.reset();
        }
        catch (ConfigurationOverlapException e) {
            overlapException = e;
        }
        neighborsUpdated();
        if (overlapException != null) {
            throw overlapException;
        }
        if (Debug.ON && Debug.DEBUG_NOW && Debug.thisPhase(phase)) {
            debugPair = Debug.getAtoms(phase);
        }
    }

    /**
     * Do an upList call for each atom and reconstruct the event list.
     */
    public void neighborsUpdated() {
        if(!initialized) return;
        super.neighborsUpdated();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom atom = atomIterator.nextAtom();
            agents[atom.getGlobalIndex()].resetCollision();
        }
        upList.setTargetAtom(null);
        collisionHandlerUp.reset();
        collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
        potential.calculate(phase, upList, collisionHandlerUp); //assumes only one phase
        eventList.reset();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Agent agent = agents[atomIterator.nextAtom().getGlobalIndex()];
            if (agent.collisionPotential != null) {
                eventList.add(agent.eventLinker);
            }
        }
    }

    /**
     * Updates collision times appropriately after scaling momenta.
     */
    protected double scaleMomenta() {
        double s = super.scaleMomenta();
        double rs = 1.0/s;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            agents[atomIterator.nextAtom().getGlobalIndex()].eventLinker.sortKey *= rs;
        }
        // don't need to update eventTree because event order didn't change
        return s;
    }

    /**
     * Updates collision times appropriately after randomizing momenta
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomenta() {
        super.randomizeMomenta();
        // super.randomizeMomenta recalculates the kinetic energy and doesn't
        // change the potential energy, so just act like the neighbor lists were
        // updated
        neighborsUpdated();
    }
    
    /**
     * Updates collision time appropriately after randomizing momentum
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomentum(Atom atom) {
        super.randomizeMomentum(atom);
        updateAtom(atom);
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
    public void setNullPotential(Potential1 nullPotential) {
        this.nullPotential = (PotentialHard)nullPotential;
        if (nullPotential != null) {
            nullPotential.setPhase(phase);
        }
    }
    
    //collision handler is passed to the potential and is notified of each collision
    //the potential detects.  The collision handler contains the necessary logic to
    //process this information so that the collision lists are kept up to date.
    //The potential should call the handler's setPotential method, with itself as the
    //argument, before beginning to detect collisions. 

    //the up-handler has the logic of the Allen & Tildesley upList subroutine
    //sets collision time of given atom to minimum value for collisions with all atoms uplist of it
    private static final class CollisionHandlerUp extends PotentialCalculation {
        double minCollisionTime;
        IntegratorHard.Agent aia;
        Atom atom1;
        double collisionTimeStep;
        private Agent[] integratorAgents;

        /**
         * resets the "atom" held by this class, ensuring the method will
         * work even if (coincidentally) the atom is the same as the same 
         * one as the last time the method was called.
         */
        public void reset() {
            atom1 = null;
        }
        
        public void setAgents(Agent[] newAgents) {
            integratorAgents = newAgents;
        }

        /**
         * sets the atom whose collision time is to be calculated
         */
        public void setAtom(Atom a) {
            atom1 = a;
            aia = integratorAgents[a.getGlobalIndex()];
            minCollisionTime = aia.collisionTime();
        }//end of setAtom

        //atom pair
        public void doCalculation(AtomsetIterator iterator, Potential potential) {
            final boolean notPairIterator = (iterator.nBody() != 2);
            iterator.reset();
            PotentialHard pHard = (PotentialHard)potential;
            while (iterator.hasNext()) {
                AtomSet atoms = iterator.next();
                if(atoms.getAtom(0) != atom1) setAtom(atoms.getAtom(0)); //need this if doing minimum collision time calculation for more than one atom
                double collisionTime = pHard.collisionTime(atoms,collisionTimeStep);
                if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atoms)) || Debug.allAtoms(atoms))) {
                    System.out.println("collision up time "+collisionTime+" for atom "+atoms+" "+pHard.getClass());
                }
                if(collisionTime < minCollisionTime) {
                    if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atoms))) {
                        System.out.println("setting up time "+collisionTime+" for atom "+atoms);
                    }
                    minCollisionTime = collisionTime;
                    aia.setCollision(collisionTime, notPairIterator ? null : atoms.getAtom(1), pHard);
                }//end if
            }
        }//end of calculate(AtomPair...

    } //end of collisionHandlerUp

	//the down-handler has the logic of the Allen & Tildesley downList subroutine
	//sets collision times of atoms downlist of given atom to minimum of their current
	//value and their value with given atom
	private static final class CollisionHandlerDown extends PotentialCalculation {
        double collisionTimeStep;
        final TreeList eventList;
        private Agent[] integratorAgents;
        CollisionHandlerDown(TreeList list) {
            eventList = list;
        }
        
        public void setAgents(Agent[] newAgents) {
            integratorAgents = newAgents;
        }
        
		public void doCalculation(AtomsetIterator iterator, Potential potential) {
			if (potential.nBody() != 2) return;
			iterator.reset();
            PotentialHard pHard = (PotentialHard)potential;
			while (iterator.hasNext()) {
				AtomPair atomPair = (AtomPair)iterator.next();
				double collisionTime = pHard.collisionTime(atomPair,collisionTimeStep);
                if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atomPair)))) {
                    System.out.println("collision down time "+collisionTime+" for atoms "+atomPair+" "+pHard.getClass());
                }
				if(collisionTime < Double.POSITIVE_INFINITY) {
					Agent aia = integratorAgents[atomPair.atom0.getGlobalIndex()];
					if(collisionTime < aia.collisionTime()) {
						if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atomPair))) {
							System.out.println("setting down time "+collisionTime+" for atoms "+atomPair);
						}
                        if (aia.collisionPotential != null) {
                            aia.eventLinker.remove();
                        }
                        aia.setCollision(collisionTime, atomPair.atom1, pHard);
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
    private static final class ReverseCollisionHandler extends PotentialCalculation {
        final AtomArrayList listToUpdate;
        private Agent[] integratorAgents;
        
        ReverseCollisionHandler(AtomArrayList list) {
            listToUpdate = list;
        }
        
        public void setAgents(Agent[] newAgents) {
            integratorAgents = newAgents;
        }
        
        public void doCalculation(AtomsetIterator iterator, Potential p) {
            if (iterator.nBody() != 2) return;
            iterator.reset();
            // look for pairs in which pair[0] is the collision partner of pair[1]
            while (iterator.hasNext()) {
                AtomPair pair = (AtomPair)iterator.next();
                Atom aPartner = integratorAgents[pair.atom0.getGlobalIndex()].collisionPartner();
                if (Debug.ON && Debug.DEBUG_NOW && ((Debug.allAtoms(pair) && Debug.LEVEL > 1) || (Debug.anyAtom(pair) && Debug.LEVEL > 2))) {
                    System.out.println(pair.atom1+" thought it would collide with "+aPartner);
                }
                if(aPartner == pair.atom1) {
                    if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(pair) || Debug.LEVEL > 2)) {
                        System.out.println("Will update "+pair.atom0+" because it wanted to collide with "+aPartner);
                    }
                    listToUpdate.add(pair.atom0);
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
    private static class CollisionListenerLinker implements java.io.Serializable {
        CollisionListenerLinker next;
        final CollisionListener listener;
        CollisionListenerLinker(CollisionListener listener, CollisionListenerLinker next) {
            this.listener = listener;
            this.next = next;
        }
    }//end of CollisionListenerLinker

    public Class getAgentClass() {
        return Agent.class;
    }
    
    /**
	 * Produces the Agent for this integrator. This method gets
     * called by the agentManager, which allocates/deallocates 
     * agents as needed.
	 */
    public Object makeAgent(Atom a) {
        if (!a.type.isLeaf()) {
            return null;
        }
        return new Agent(a,this);
    }
    
    // don't need to remove the agent from the event list because reset will
    // get called and that will totally clear the event list
    public void releaseAgent(Object agent, Atom atom) {}
 
    /**
     * Agent defined by this integrator.
     * Holds information about the time to next collision (considering only atoms
     * uplist of the atom), the collisionPartner (atom that this one will be colliding with),
     * and a handle to the potential governing the interactions of this atom and its collision partner
     * (this is kept for convenience, so it won't have to be determined again when the collision
     * is processed).
     */
    //Do not use encapsulation since the fields are for referencing by the integrator
    public static class Agent implements Serializable {  //need public so to use with instanceof
        public Atom atom, collisionPartner;
        public PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        public TreeLinker eventLinker;
        private final IntegratorHard integrator;
        
        public Agent(Atom a, IntegratorHard integrator) {
            this.integrator = integrator;
            atom = a;
            eventLinker = new TreeLinker(this);
            eventLinker.sortKey = Double.POSITIVE_INFINITY;
        }
        
        public String toString() {
            return "Collider: "+atom+"; Partner: "+collisionPartner+"; Potential: "+collisionPotential;
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
            collisionPotential = integrator.nullPotential;
            if (collisionPotential != null) {
                eventLinker.sortKey = collisionPotential.collisionTime(atom,integrator.collisionTimeStep);
            }
            else {
                eventLinker.sortKey = Double.POSITIVE_INFINITY;
            }
            collisionPartner = null;
        }
        
        /**
         * Sets parameters associated with next two-body collision of this atom with another atom.
         *
         * @param time    time to collision of this agent's atom with an atom uplist of it
         * @param partner the atom this one will collide with next
         * @param p       the potential for interactions between this atom and its collision partner
         */
        public final void setCollision(double time, Atom partner, PotentialHard p) {
            collisionPartner = partner;
            collisionPotential = p;
            eventLinker.sortKey = time;
        }


        /**
         * Decreases the recorded time to collision of this atom
         * This action is performed when the atom is advanced without a collision
         */
        public final void decrementCollisionTime(double interval) {
            eventLinker.sortKey -= interval;
        }
        /**
         * Accessor method for the time to next collision of this atom
         */
        public final double collisionTime() {return eventLinker.sortKey;}
    }//end of Agent
    
    public interface CollisionListener {
        public void collisionAction(Agent colliderAgent);
    }

}

