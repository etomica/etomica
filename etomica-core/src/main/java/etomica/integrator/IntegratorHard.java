/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.api.*;
import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.exception.ConfigurationOverlapException;
import etomica.nbr.PotentialMasterHybrid;
import etomica.nbr.list.INeighborListListener;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialHard;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.util.Debug;
import etomica.util.TreeLinker;
import etomica.util.TreeList;

/**
 * Integrator for hard potentials.
 * Integrates equations of motion through time by advancing atoms from one
 * collision to the next.  Determination of time of collision and
 * implementation of collision dynamics is handled by the potential between the
 * atoms and stored in an Agent.
 *
 * @author David Kofke
 *
 */
public class IntegratorHard extends IntegratorMD
    implements INeighborListListener, AgentSource<IntegratorHard.Agent>, AtomTypeAgentManager.AgentSource {

    private static final long serialVersionUID = 1L;
    //handle to the integrator agent holding information about the next collision
    protected IntegratorHard.Agent colliderAgent;
    //first of a linked list of objects (typically meters) that are called each time a collision is processed
    protected CollisionListenerLinker collisionListenerHead = null;
    private final AtomPair pair;
    private final AtomSetSinglet singlet;
    private double minDelta;
    private AtomPair debugPair;

    protected final IteratorDirective upList = new IteratorDirective(IteratorDirective.Direction.UP);
    protected final IteratorDirective downList = new IteratorDirective(IteratorDirective.Direction.DOWN);
    protected final AtomArrayList listToUpdate = new AtomArrayList();
    protected final TreeList eventList = new TreeList();

    protected final ReverseCollisionHandler reverseCollisionHandler;
    protected final CollisionHandlerUp collisionHandlerUp;
    protected final CollisionHandlerDown collisionHandlerDown;

    protected double collisionTimeStep;
    protected long collisionCount;

    protected AtomLeafAgentManager<Agent> agentManager;
    protected final AtomTypeAgentManager nullPotentialManager;
    protected int handlingEvent;

    public IntegratorHard(Simulation sim, PotentialMaster potentialMaster, Space _space) {
        this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }

    public IntegratorHard(Simulation sim, PotentialMaster potentialMaster, IRandom random,
                          double timeStep, double temperature, Space _space) {
        super(potentialMaster,random,timeStep,temperature, _space);
        pair = new AtomPair();
        singlet = new AtomSetSinglet();
        reverseCollisionHandler = new ReverseCollisionHandler(listToUpdate);
        reverseCollisionHandler.setAgentManager(agentManager);
        collisionHandlerUp = new CollisionHandlerUp();
        collisionHandlerUp.setAgentManager(agentManager);
        collisionHandlerDown = new CollisionHandlerDown(eventList);
        collisionHandlerDown.setAgentManager(agentManager);
        if (sim != null) {
            nullPotentialManager = new AtomTypeAgentManager(this);
            nullPotentialManager.init(sim);
        }
        else {
            nullPotentialManager = null;
        }

    }
    
    public void setBox(Box newBox) {
        if (box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
            if(this.potentialMaster instanceof PotentialMasterList) {
                ((PotentialMasterList)this.potentialMaster).getNeighborManager(box).getEventManager().removeListener(this);
            }
        }
        super.setBox(newBox);
        agentManager = new AtomLeafAgentManager<Agent>(this,newBox,Agent.class);
        collisionHandlerUp.setAgentManager(agentManager);
        collisionHandlerDown.setAgentManager(agentManager);
        reverseCollisionHandler.setAgentManager(agentManager);

        if (nullPotentialManager != null) {
            AtomTypeAgentManager.AgentIterator iterator = nullPotentialManager.makeIterator();
            iterator.reset();
            while (iterator.hasNext()) {
                ((PotentialHard)iterator.next()).setBox(newBox);
            }
        }
        if(this.potentialMaster instanceof PotentialMasterList) {
            ((PotentialMasterList)this.potentialMaster).getNeighborManager(box).getEventManager().addListener(this);
        }
        else if (this.potentialMaster instanceof PotentialMasterHybrid) {
            ((PotentialMasterHybrid)this.potentialMaster).getNeighborManager(box).getEventManager().addListener(this);
        }
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
    public void doStepInternal() {
        if (Double.isInfinite(currentPotentialEnergy)) {
            // we were overlapped at some point.  try recalculating the PE now
            // so we can start re-tracking the PE once we aren't overlapped.
            currentPotentialEnergy = meterPE.getDataAsScalar();
        }
        super.doStepInternal();
        findNextCollider();
        collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.POSITIVE_INFINITY;
        double oldTime = 0;
        while(collisionTimeStep < timeStep) {//advance to collision if occurs before remaining interval
            IAtomList atoms;
            if (colliderAgent.collisionPartner() != null) {
                atoms = pair;
                pair.atom0 = colliderAgent.atom();
                pair.atom1 = colliderAgent.collisionPartner();
            }
            else {
                atoms = singlet;
                singlet.atom = colliderAgent.atom();
            }
            if (collisionTimeStep - oldTime < minDelta) {
                System.out.println("diff "+(collisionTimeStep - oldTime)+" minDelta "+minDelta);
                System.out.println("previous collision occured after current one");
                System.out.println("previous time: "+oldTime+" current time: "+collisionTimeStep);
                System.out.println("collision for "+atoms+" potential "+colliderAgent.collisionPotential.getClass());
                if (atoms instanceof AtomPair) {
                    Vector dr = space.makeVector();
                    Vector dv = space.makeVector();

                    IAtomKinetic atom0 = (IAtomKinetic)pair.atom0;
                    IAtomKinetic atom1 = (IAtomKinetic)pair.atom1;
                    dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
                    
                    dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
                    box.getBoundary().nearestImage(dr);

                    dr.PEa1Tv1(oldTime,dv);
                    System.out.println("distance at last collision time was "+Math.sqrt(dr.squared()));
                    dr.PEa1Tv1(collisionTimeStep-oldTime,dv);
                    System.out.println("distance now "+Math.sqrt(dr.squared()));
                }
                throw new RuntimeException("this simulation is not a time machine");
            }
            if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 0 &&((Debug.LEVEL > 1 && Debug.thisBox(box)) || Debug.anyAtom(atoms))) {
                System.out.println("collision between atoms "+atoms+" at "+collisionTimeStep+" with "+colliderAgent.collisionPotential.getClass());
            }
            if (Debug.ON && Debug.DEBUG_NOW && Debug.thisBox(box) && (debugPair = Debug.getAtoms(box)) != null) {
                if (Debug.LEVEL > 1) {
                    IAtom a = debugPair.atom0;
                    if (a != null) {
                        Agent agent = agentManager.getAgent(a);
                        System.out.println(a+" collision with "+agent.collisionPartner+" "+agent.collisionPotential+" at "+agent.collisionTime());
                    }
                    a = debugPair.atom1;
                    if (a != null) {
                        Agent agent = agentManager.getAgent(a);
                        System.out.println(a+" collision with "+agent.collisionPartner+" "+agent.collisionPotential+" at "+agent.collisionTime());
                    }
                }
                if (debugPair.atom0 != null && debugPair.atom1 != null && !(debugPair.atom0 instanceof IMolecule && debugPair.atom1 instanceof IMolecule)) {
                    Vector dr = space.makeVector();
                    Vector dv = space.makeVector();

                    IAtomKinetic atom0 = (IAtomKinetic)debugPair.atom0;
                    IAtomKinetic atom1 = (IAtomKinetic)debugPair.atom1;
                    dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
                    
                    dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());

                    dr.PEa1Tv1(collisionTimeStep,dv);
                    box.getBoundary().nearestImage(dr);
                    double r2 = dr.squared();
                    if (Debug.LEVEL > 1 || Math.sqrt(r2) < Debug.ATOM_SIZE-1.e-11) {
                        System.out.println("distance between "+debugPair+" is "+Math.sqrt(r2));
                        if (Debug.LEVEL > 2 || Math.sqrt(r2) < Debug.ATOM_SIZE-1.e-11) {
                            dr.Ea1Tv1(collisionTimeStep,atom0.getVelocity());
                            dr.PE(atom0.getPosition());
                            System.out.println(atom0+" coordinates "+dr);
                            dr.Ea1Tv1(collisionTimeStep,atom1.getVelocity());
                            dr.PE(atom1.getPosition());
                            System.out.println(atom1+" coordinates "+dr);
                        }
                    }
                    else if (Debug.LEVEL > 2 && !(debugPair.atom0 instanceof IMolecule)) {
                        dr.Ea1Tv1(collisionTimeStep,((IAtomKinetic)debugPair.atom0).getVelocity());
                        dr.PE(debugPair.atom0.getPosition());
                        System.out.println(debugPair.atom0+" coordinates "+dr);
                    }
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
            if (atoms.getAtomCount() == 1) {
                updateAtom(atoms.getAtom(0));
            }
            else {
                updateAtoms((AtomPair)atoms);
            }
            findNextCollider(); //this sets colliderAgent for the next collision
            if (Debug.ON && colliderAgent != null && colliderAgent.atom == atoms.getAtom(0) && (atoms.getAtomCount() == 2 && colliderAgent.collisionPartner == atoms.getAtom(1))
                    && colliderAgent.collisionTime() == collisionTimeStep) {
                throw new RuntimeException("repeating collision "+atoms+" "+collisionTimeStep);
            }
            
            oldTime = collisionTimeStep;
            collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.POSITIVE_INFINITY;
        } 

        advanceAcrossTimeStep(timeStep);
        collisionTimeStep = 0.0;
        if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.thisBox(box)) {
            eventList.check();
            double PE = meterPE.getDataAsScalar();
            if (Math.abs((PE - currentPotentialEnergy)/(PE+currentPotentialEnergy)) > 1.e-9
                    && Math.abs(PE - currentPotentialEnergy) > 1.e-9) {
                throw new RuntimeException("potential energy of "+box+" is wrong. it's actually "+PE+" but I thought it was "+currentPotentialEnergy);
            }
            double KE = meterKE.getDataAsScalar();
            if (Math.abs(KE - currentKineticEnergy) > 1.e-8) {
                throw new RuntimeException("kinetic energy of "+box+" is wrong. it's actually "+KE+" but I thought it was "+currentKineticEnergy);
            }
        }

        if(isothermal) doThermostatInternal();
    }//end of doStep

    public long getCollisionCount() {
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
        potentialMaster.calculate(box, downList, reverseCollisionHandler);
        downList.setTargetAtom(colliders.atom1);
        potentialMaster.calculate(box, downList, reverseCollisionHandler);

        // this would update collider0 as well since it wanted to collide with
        // atom1.  But we to full reset it, so remove it from the (hopefully
        // small) list.
        listToUpdate.remove(listToUpdate.indexOf(colliders.atom0));
        processReverseList();

        Agent agent = agentManager.getAgent(colliders.atom0);
        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollisionFull();
        upList.setTargetAtom(colliders.atom0);
        collisionHandlerUp.setAtom(colliders.atom0);
        collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
        potentialMaster.calculate(box, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        downList.setTargetAtom(colliders.atom0);
        collisionHandlerDown.collisionTimeStep = this.collisionTimeStep;
        potentialMaster.calculate(box, downList, collisionHandlerDown);

        agent = agentManager.getAgent(colliders.atom1);
        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollisionFull();
        upList.setTargetAtom(colliders.atom1);
        collisionHandlerUp.setAtom(colliders.atom1);
        collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
        potentialMaster.calculate(box, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        downList.setTargetAtom(colliders.atom1);
        collisionHandlerDown.collisionTimeStep = this.collisionTimeStep;
        potentialMaster.calculate(box, downList, collisionHandlerDown);
    }

    /**
     * updates collision time for a single atom (and any atom that might a 
     * collision with the atom)
     */
    protected void updateAtom(IAtom a) {
        Agent agent = agentManager.getAgent(a);
        if(agent == null) return;

        listToUpdate.clear();

        downList.setTargetAtom(a);
        potentialMaster.calculate(box, downList, reverseCollisionHandler);
        processReverseList();

        if (agent.collisionPotential != null) {
            agent.eventLinker.remove();
        }
        agent.resetCollisionFull();
        upList.setTargetAtom(a);
        collisionHandlerUp.setAtom(a);
        collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
        potentialMaster.calculate(box, upList, collisionHandlerUp);
        if (agent.collisionPotential != null) {
            eventList.add(agent.eventLinker);
        }
        collisionHandlerDown.collisionTimeStep = this.collisionTimeStep;
        potentialMaster.calculate(box, downList, collisionHandlerDown);
    }

    /**
     * Recalculate collision times (up) for any atoms atoms in the "reverse"
     * list of atoms that thought they would collide with an atom that's
     * changed.
     */
    protected void processReverseList() {
        int size = listToUpdate.getAtomCount();
        for (int i=0; i<size; i++) {
            IAtom reverseAtom = listToUpdate.getAtom(i);
            Agent agent = agentManager.getAgent(reverseAtom);
            if (agent.collisionPotential != null) {
                agent.eventLinker.remove();
            }
            // reset collision, but not a "full" reset
            // this atom thought it would collide with something and now it
            // won't.  We should rever to the "null" collision time it had
            // the last time it had a real collision of its own.
            agent.resetCollision();
            upList.setTargetAtom(reverseAtom);
            collisionHandlerUp.collisionTimeStep = this.collisionTimeStep;
            collisionHandlerUp.setAtom(reverseAtom);
            potentialMaster.calculate(box, upList, collisionHandlerUp);
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
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            agentManager.getAgent(a).decrementCollisionTime(tStep);
			a.getPosition().PEa1Tv1(tStep,a.getVelocity());
		}
	}

    public void reset() {
        ConfigurationOverlapException overlapException = null;
        try {
            super.reset();
        }
        catch (ConfigurationOverlapException e) {
            overlapException = e;
        }
        colliderAgent = null;
        resetCollisionTimes();
        if (overlapException != null) {
            throw overlapException;
        }
        if (Debug.ON && Debug.DEBUG_NOW && Debug.thisBox(box)) {
            debugPair = Debug.getAtoms(box);
        }
    }

    public void neighborListNeighborsUpdated() {
        handlingEvent++;
        resetCollisionTimes();
        handlingEvent--;
    }
    
    public void boxMoleculeAdded(IBoxMoleculeEvent e) {
        handlingEvent++;
        super.boxMoleculeAdded(e);
        handlingEvent--;
    }

    /**
     * Do an upList call for each atom and reconstruct the event list.
     */
    public void resetCollisionTimes() {
        if(!initialized) return;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom atom = leafList.getAtom(iLeaf);
            agentManager.getAgent(atom).resetCollisionFull();
        }
        upList.setTargetAtom(null);
        collisionHandlerUp.reset();
        collisionHandlerUp.collisionTimeStep = 0;
        potentialMaster.calculate(box, upList, collisionHandlerUp); //assumes only one box
        eventList.reset();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom atom = leafList.getAtom(iLeaf);
            Agent agent = agentManager.getAgent(atom);
            if (agent.collisionPotential != null) {
                eventList.add(agent.eventLinker);
            }
        }
    }

    /**
     * Updates collision times appropriately after scaling momenta.
     */
    protected void scaleMomenta() {
        super.scaleMomenta();
        // super.scaleMomenta alters the velocities, so we need to 
        // recalculate collision times
        resetCollisionTimes();
    }

    /**
     * Updates collision times appropriately after randomizing momenta
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomenta() {
        super.randomizeMomenta();
        // super.randomizeMomenta alters the velocities, so we need to 
        // recalculate collision times
        resetCollisionTimes();
    }

    /**
     * Updates collision times appropriately after randomizing momenta
     * as part of the Andersen scaling thermostat.
     */
    protected void randomizeTotalKE() {
        super.randomizeTotalKE();
        // super.randomizeMomenta alters the velocities, so we need to 
        // recalculate collision times
        resetCollisionTimes();
    }

    /**
     * Updates collision time appropriately after randomizing momentum
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomentum(IAtomKinetic atom) {
        super.randomizeMomentum(atom);
        if (handlingEvent == 0) {
            updateAtom(atom);
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
     * @return Returns the nullPotential.
     */
    public PotentialHard getNullPotential(IAtomType atomType) {
        return (PotentialHard)nullPotentialManager.getAgent(atomType);
    }

    /**
     * Sets the "null" potential for the given AtomType.  The null potential is
     * a 1-body potential used by IntegratorHard to handle Atoms that wrap
     * around periodic boundaries when neighbor listing is not used.
     */
    public void setNullPotential(PotentialHard nullPotential, IAtomType type) {
        // if nullPotentialManager is null, it's because you passed a null
        // ISimulation when you constructed this class
        nullPotentialManager.setAgent(type, nullPotential);
        if (nullPotential != null && box != null) {
            nullPotential.setBox(box);
        }
        if (box != null) {
            // inefficient -- we could loop over molecules of type.getSpecies()
            // and then their children.  oh well.
            IAtomList leafList = box.getLeafList();
            int nLeaf = leafList.getAtomCount();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom atom = leafList.getAtom(iLeaf);
                if (atom.getType() == type) {
                    agentManager.getAgent(atom).setNullPotential(nullPotential);
                }
            }
        }
    }

    //collision handler is passed to the potential and is notified of each collision
    //the potential detects.  The collision handler contains the necessary logic to
    //process this information so that the collision lists are kept up to date.
    //The potential should call the handler's setPotential method, with itself as the
    //argument, before beginning to detect collisions. 

    //the up-handler has the logic of the Allen & Tildesley upList subroutine
    //sets collision time of given atom to minimum value for collisions with all atoms uplist of it
    protected static final class CollisionHandlerUp implements PotentialCalculation {
        double minCollisionTime;
        IntegratorHard.Agent aia;
        IAtom atom1;
        double collisionTimeStep;
        private AtomLeafAgentManager<Agent> integratorAgentManager;

        /**
         * resets the "atom" held by this class, ensuring the method will
         * work even if (coincidentally) the atom is the same as the same 
         * one as the last time the method was called.
         */
        public void reset() {
            atom1 = null;
        }

        public void setAgentManager(AtomLeafAgentManager<Agent> newAgentManager) {
            integratorAgentManager = newAgentManager;
        }

        /**
         * sets the atom whose collision time is to be calculated
         */
        public void setAtom(IAtom a) {
            atom1 = a;
            aia = integratorAgentManager.getAgent(a);
            minCollisionTime = aia.collisionTime();
        }//end of setAtom

        //atom pair
        public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
            PotentialHard pHard = (PotentialHard)potential;
            if(atoms.getAtom(0) != atom1) setAtom(atoms.getAtom(0)); //need this if doing minimum collision time calculation for more than one atom
            double collisionTime = pHard.collisionTime(atoms,collisionTimeStep);
            if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atoms)) || Debug.allAtoms(atoms))) {
                System.out.println("collision up time "+collisionTime+" for atom "+atoms+" "+pHard.getClass());
            }
            if(collisionTime < minCollisionTime) {
                if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && (Debug.LEVEL > 2 || Debug.anyAtom(atoms))) {
                    System.out.println("setting up time "+collisionTime+" for atom "+atoms);
                }
                minCollisionTime = collisionTime;
                aia.setCollision(collisionTime, (atoms.getAtomCount() == 2) ? atoms.getAtom(1) : null, pHard);
            }//end if
        }//end of calculate(AtomPair...

    } //end of collisionHandlerUp

	//the down-handler has the logic of the Allen & Tildesley downList subroutine
	//sets collision times of atoms downlist of given atom to minimum of their current
	//value and their value with given atom
	private static final class CollisionHandlerDown implements PotentialCalculation, java.io.Serializable {
        private static final long serialVersionUID = 1L;
        double collisionTimeStep;
        final TreeList eventList;
        private AtomLeafAgentManager<Agent> integratorAgentManager;
        CollisionHandlerDown(TreeList list) {
            eventList = list;
        }

        public void setAgentManager(AtomLeafAgentManager<Agent> newAgentManager) {
            integratorAgentManager = newAgentManager;
        }

		public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
			if (atoms.getAtomCount() != 2) return;
            PotentialHard pHard = (PotentialHard)potential;

			double collisionTime = pHard.collisionTime(atoms,collisionTimeStep);
            if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atoms)))) {
                System.out.println("collision down time "+collisionTime+" for atoms "+atoms+" "+pHard.getClass());
            }
			if(collisionTime < Double.POSITIVE_INFINITY) {
				Agent aia = integratorAgentManager.getAgent(atoms.getAtom(0));
				if(collisionTime < aia.collisionTime()) {
					if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atoms))) {
						System.out.println("setting down time "+collisionTime+" for atoms "+atoms);
					}
                    if (aia.collisionPotential != null) {
                        aia.eventLinker.remove();
                    }
                    aia.setCollision(collisionTime, atoms.getAtom(1), pHard);
                    eventList.add(aia.eventLinker);
				}
			}
		}
	}

    /**
     * A PotentialCalculation to find atoms that thought they would collide
     * with an atom.  The iterator should return an atom and its "down"
     * neighbors.  
     */
    private static final class ReverseCollisionHandler implements PotentialCalculation, java.io.Serializable {
        private static final long serialVersionUID = 1L;
        final AtomArrayList listToUpdate;
        private AtomLeafAgentManager<Agent> integratorAgentManager;
        
        ReverseCollisionHandler(AtomArrayList list) {
            listToUpdate = list;
        }
        
        public void setAgentManager(AtomLeafAgentManager<Agent> newAgentManager) {
            integratorAgentManager = newAgentManager;
        }
        
        public void doCalculation(IAtomList pair, IPotentialAtomic p) {
            if (pair.getAtomCount() != 2) return;
            // look for pairs in which pair[0] is the collision partner of pair[1]
            IAtom aPartner = integratorAgentManager.getAgent(pair.getAtom(0)).collisionPartner();
            if (Debug.ON && Debug.DEBUG_NOW && ((Debug.allAtoms(pair) && Debug.LEVEL > 1) || (Debug.anyAtom(pair) && Debug.LEVEL > 2))) {
                System.out.println(pair.getAtom(0)+" thought it would collide with "+aPartner);
            }
            if(aPartner == pair.getAtom(1)) {
                if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(pair) || Debug.LEVEL > 2)) {
                    System.out.println("Will update "+pair.getAtom(0)+" because it wanted to collide with "+aPartner);
                }
                listToUpdate.add(pair.getAtom(0));
            }
        }
    }

    /**
     * Class used to construct a linked list of collision listeners
     */
    private static class CollisionListenerLinker implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        CollisionListenerLinker next;
        final CollisionListener listener;
        CollisionListenerLinker(CollisionListener listener, CollisionListenerLinker next) {
            this.listener = listener;
            this.next = next;
        }
    }//end of CollisionListenerLinker

    /**
	 * Produces the Agent for this integrator. This method gets
     * called by the agentManager, which allocates/deallocates 
     * agents as needed.
	 */
    public Agent makeAgent(IAtom a, Box agentBox) {
        Agent agent = new Agent(a, this);
        if (nullPotentialManager != null) {
            agent.setNullPotential((PotentialHard)nullPotentialManager.getAgent(a.getType()));
        }
        return agent;
    }

    // don't need to remove the agent from the event list because reset will
    // get called and that will totally clear the event list
    public void releaseAgent(Agent agent, IAtom atom, Box agentBox) {}

    public Class getSpeciesAgentClass() {
        return PotentialHard.class;
    }

    public Object makeAgent(IAtomType type) {
        return null;
    }

    public void releaseAgent(Object agent, IAtomType type) {
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
    public static class Agent {  //need public so to use with instanceof
        public IAtom atom, collisionPartner;
        public PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        public TreeLinker eventLinker;
        protected PotentialHard nullPotential;
        protected AtomSetSinglet atomSetSinglet;
        protected double nullCollisionTime;
        protected final IntegratorHard integrator;

        public Agent(IAtom a, IntegratorHard integrator) {
            atom = a;
            eventLinker = new TreeLinker(this);
            eventLinker.sortKey = Double.POSITIVE_INFINITY;
            nullCollisionTime = Double.POSITIVE_INFINITY;
            this.integrator = integrator;
        }

        public void setNullPotential(PotentialHard newNullPotential) {
            nullPotential = newNullPotential;
        }

        public String toString() {
            return "Collider: "+atom+"; Partner: "+collisionPartner+"; Potential: "+collisionPotential;
        }

        public final IAtom atom() {return atom;}
        public final IAtom collisionPartner() {return collisionPartner;}

        /**
         * resets collision potential and partner.  If a null potnetial is in
         * use, the time is reset to the null potential collision time.
         * caller should remove eventLinker from the tree if needed before
         * calling this method.
         */
        public void resetCollision() {
            // events with a (non-null) potential must be in the tree
            // events in the tree must have a potential
            collisionPotential = nullPotential;
            if (collisionPotential != null) {
                eventLinker.sortKey = nullCollisionTime;
                if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.anyAtom(atomSetSinglet)) {
                    System.out.println("resetting collision time back to null collision time for "+atom+" to "+eventLinker.sortKey);
                }
            }
            else {
                eventLinker.sortKey = Double.POSITIVE_INFINITY;
            }
            collisionPartner = null;
        }

        /**
         * resets time, potential and partner.  caller should remove 
         * eventLinker from the tree if needed before calling this method.
         */
        public void resetCollisionFull() {
            // events with a (non-null) potential must be in the tree
            // events in the tree must have a potential
            collisionPotential = nullPotential;
            if (collisionPotential != null) {
                if (atomSetSinglet == null) {
                    atomSetSinglet = new AtomSetSinglet();
                }
                atomSetSinglet.atom = atom;
                eventLinker.sortKey = collisionPotential.collisionTime(atomSetSinglet,integrator.collisionTimeStep);
                nullCollisionTime = eventLinker.sortKey;
                if (Debug.ON && Debug.DEBUG_NOW && Debug.LEVEL > 1 && Debug.anyAtom(atomSetSinglet)) {
                    System.out.println("initializing null collision time for "+atom+" to "+eventLinker.sortKey);
                }
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
        public final void setCollision(double time, IAtom partner, PotentialHard p) {
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
            nullCollisionTime -= interval;
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

