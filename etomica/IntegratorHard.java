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
 
public class IntegratorHard extends IntegratorMD {

    //handle to the integrator agent holding information about the next collision
    protected IntegratorHard.Agent colliderAgent;
    //iterators for looping through atoms
    protected AtomIteratorListSimple atomIterator = new AtomIteratorListSimple();
    //first of a linked list of objects (typically meters) that are called each time a collision is processed
    protected CollisionListenerLinker collisionListenerHead = null;
    //time elapsed since reaching last timestep increment
    private double timeIncrement = 0.0;
    private Atom[] atoms;
    protected final MeterTemperature meterTemperature = new MeterTemperature();
    Space.Vector c3;
    Space.CoordinatePair cPairDebug;

	protected final IteratorDirective upList = new IteratorDirective(IteratorDirective.UP);
	protected final IteratorDirective downList = new IteratorDirective(IteratorDirective.DOWN);
	protected final CollisionHandlerUp collisionHandlerUp = new CollisionHandlerUp();
	protected final CollisionHandlerDown collisionHandlerDown = new CollisionHandlerDown();

	private Atom[] targetAtom = new Atom[1];
	
    
    public IntegratorHard(PotentialMaster potentialMaster) {
        super(potentialMaster);
        Agent.nullPotential = null; //(PotentialHard)Potential.NullPotential(sim);
        atoms = new Atom[2];
    }//end of constructor
    
    public boolean addPhase(Phase phase) {
        if(!super.addPhase(phase)) return false;
        atomIterator.setList(phase.speciesMaster.atomList);
        meterTemperature.setPhase(this.phase);
        return true;
    }
    
    public IntegratorHard.Agent colliderAgent() {
        return colliderAgent;
    }
    
    /** 
     * Steps all atoms across time interval timeStep, handling all intervening collisions.
     */
    public void doStep() {
        double collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.MAX_VALUE;
        double interval = timeStep;
        int count = 10000;
        while(collisionTimeStep < interval) {//advance to collision if occurs before remaining interval
			atoms[0] = colliderAgent.atom();
			atoms[1] = colliderAgent.collisionPartner();
        	if (collisionTimeStep < 0.0) {
        		System.out.println("previous collision occured before current one");
        		System.out.println("previous time: "+(timeStep-interval)+"current time: "+(timeStep-interval+collisionTimeStep));
        		System.out.println("collision between "+atoms[0]+" and "+atoms[1]);
        		cPairDebug = Simulation.getDefault().space.makeCoordinatePair();
        		cPairDebug.setBoundary(firstPhase.boundary());
        		cPairDebug.reset(atoms[0].coord,atoms[1].coord);
        		System.out.println("distance at last collision time was "+cPairDebug.r2());
                advanceAcrossTimeStep(collisionTimeStep);//if needing more flexibility, make this a separate method-- advanceToCollision(collisionTimeStep)
        		cPairDebug.reset();
        		System.out.println("distance now "+cPairDebug.r2());
        		throw new RuntimeException("this simulation is not a time machine");
        	}
        	if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 1 || Debug.anyAtom(atoms))) {
        		System.out.println("collision between atoms "+atoms[0]+" and "+atoms[1]+" at "+(timeStep-interval+collisionTimeStep));
        	}
            advanceAcrossTimeStep(collisionTimeStep);//if needing more flexibility, make this a separate method-- advanceToCollision(collisionTimeStep)
    		if (Debug.ON && Debug.DEBUG_NOW && Debug.ATOM1 != null && Debug.ATOM2 != null) {
    			cPairDebug = Simulation.getDefault().space.makeCoordinatePair();
    			cPairDebug.setBoundary(firstPhase.boundary());
    			Debug.checkAtoms(cPairDebug);
    		}
			if (colliderAgent.collisionPotential != null) {
				colliderAgent.collisionPotential.bump(atoms);
			}
            for(CollisionListenerLinker cll=collisionListenerHead; cll!=null; cll=cll.next) {
                cll.listener.collisionAction(colliderAgent);
            }
            updateCollisions();
            findNextCollider(); //this sets colliderAgent for the next collision
            
            interval -= collisionTimeStep;
            collisionTimeStep = (colliderAgent != null) ? colliderAgent.collisionTime() : Double.MAX_VALUE;
            if(count-- == 0) throw new RuntimeException("Unable to advance system through all collisions");
        } 
        advanceAcrossTimeStep(interval);
        if(isothermal) {
            scaleMomenta(Math.sqrt(this.temperature/meterTemperature.getDataAsScalar(firstPhase)));
        }
        
    }//end of doStep
    
   /**
	* Loops through all atoms to identify the one with the smallest value of collisionTime
	* Collision time is obtained from the value stored in the Integrator.Agent from each atom.
	*/
	protected void findNextCollider() {
		//find next collision pair by looking for minimum collisionTime
		double minCollisionTime = Double.MAX_VALUE;
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			IntegratorHard.Agent ia = (IntegratorHard.Agent)atomIterator.nextAtom().ia;
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
//		  if (partner != null) System.out.println(collider.toString()+"\t"+partner.toString());
//		  else System.out.println(collider.toString());
            
	//   Do upList for any atoms that were scheduled to collide with atoms colliding now
	//   Assumes collider and partner haven't moved in list
		atomIterator.reset();  //first atom in first cell
		while(atomIterator.hasNext()) {
            
			//to be fixed...
			//as constructed using AtomIteratorListSimple with speciesMaster.atomList as the basis,
			// this atom iterator is not going through the atoms in sequence
			Atom a = atomIterator.nextAtom();
			if(a == collider) { //finished with atoms before collider...
				if(partner == null) break;  //and there is no partner, so we're done, or...
				continue;                   //...else just skip this atom and continue with loop
			}
			if(a == partner) break; //finished with atoms before partner; we're done
			IntegratorHard.Agent aAgent = (IntegratorHard.Agent)a.ia;
			Atom aPartner = aAgent.collisionPartner();
			if(aPartner != null && (aPartner == partner || aPartner == collider)) {//aPartner != null handles case where aPartner and partner are both null
				aAgent.resetCollision();
				targetAtom[0] = a;
				upList.setTargetAtoms(targetAtom);
				potential.calculate(firstPhase, upList, collisionHandlerUp.setAtom(a));
			}
		}//end while

		colliderAgent.resetCollision();
		targetAtom[0] = collider;
		upList.setTargetAtoms(targetAtom);
		downList.setTargetAtoms(targetAtom);
		potential.calculate(firstPhase, upList, collisionHandlerUp.setAtom(collider));
		potential.calculate(firstPhase, downList, collisionHandlerDown);
		if(partner != null) {
			((IntegratorHard.Agent)partner.ia).resetCollision();
			targetAtom[0] = partner;
			upList.setTargetAtoms(targetAtom);
			downList.setTargetAtoms(targetAtom);
			potential.calculate(firstPhase, upList, collisionHandlerUp.setAtom(partner));
			potential.calculate(firstPhase, downList, collisionHandlerDown);
		}

	}//end of processCollision
    
	/**
	 * updates collision time for a single atom (and any atom that might a collision
	 * with the atom)
	 * @param a
	 */
	protected void updateAtom(Atom a) {
		Agent agent = (Agent)a.ia;
		agent.resetCollision();
		targetAtom[0] = a;
		upList.setTargetAtoms(targetAtom);
		downList.setTargetAtoms(targetAtom);
		potential.calculate(firstPhase, upList, collisionHandlerUp.setAtom(a));
		potential.calculate(firstPhase, downList, collisionHandlerDown);
	}

	/**
	* Advances all atom coordinates by tStep, without any intervening collisions.
	* Uses free-flight kinematics.
	*/
	protected void advanceAcrossTimeStep(double tStep) {
        
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			Atom a = atomIterator.nextAtom();
	 //   for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
			((Agent)a.ia).decrementCollisionTime(tStep);
			if(a.coord.isStationary()) {continue;}  //skip if atom is stationary
			a.coord.freeFlight(tStep);
	//     a.translateBy(tStep*a.rm(),a.momentum());
		}
	}

	/**
	* Sets up the integrator for action.
	* Do an upList call for each atom and find the next collider
	*/
	protected void reset() {
		if(isothermal) scaleMomenta(Math.sqrt(this.temperature/meterTemperature.getDataAsScalar(firstPhase)));
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			((IntegratorHard.Agent)atomIterator.nextAtom().ia).resetCollision();
		}
		targetAtom[0] = null;
		upList.setTargetAtoms(targetAtom);
		potential.calculate(firstPhase, upList, collisionHandlerUp); //assumes only one phase
		findNextCollider();
	}
                
    /**
    * The simulation time elapsed since the start of the integration.
    * Cannot be reset to zero.  
    * Modified from superclass method to return correct value if in the middle
    * of processing collisions between time steps.
    */
//    public double elapsedTime() {return super.elapsedTime() + timeIncrement;}

              
    /**
    * Crude method to enforce constant-temperature constraint
    * Scales momenta of all atoms by a constant factor so that 
    * phase adheres to setpoint temperature.  Updates collision times appropriately.
    */
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            a.coord.momentum().TE(s); //scale momentum
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
        
		public CollisionHandlerUp setAtom(Atom a) {
			atom1 = a;
			aia = (IntegratorHard.Agent)a.ia;
			minCollisionTime = aia.collisionTime();
			return this;
		}//end of setAtom
        
		//atom pair
		public void doCalculation(AtomsetIterator iterator, Potential potential) {
//			System.out.println("collisionhandlerup "+pair.toString());
			iterator.reset();
			PotentialHard pHard = (PotentialHard)potential; 
			while (iterator.hasNext()) {
				Atom[] atoms = iterator.next();
				if(atoms[0] != atom1) setAtom(atoms[0]); //need this if doing minimum collision time calculation for more than one atom
				double collisionTime = pHard.collisionTime(atoms);
				if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atoms)))) {
					System.out.println("collision time "+collisionTime+" for atom "+atoms[0]+(atoms.length > 1 ? " with "+atoms[1] : ""));
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
	private static final class CollisionHandlerDown extends PotentialCalculation {
		public void doCalculation(AtomsetIterator iterator, Potential potential) {
			if (potential.nBody() != 2) return;
			iterator.reset();
			PotentialHard pHard = (PotentialHard)potential; 
			while (iterator.hasNext()) {
				Atom[] atoms = iterator.next();
				double collisionTime = pHard.collisionTime(atoms);
				if(collisionTime < Double.MAX_VALUE) {
					IntegratorHard.Agent aia = (IntegratorHard.Agent)atoms[1].ia;
					if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || (Debug.LEVEL > 1 && Debug.anyAtom(atoms)))) {
						System.out.println("collision time "+collisionTime+" for atom "+atoms[1]+(atoms.length > 1 ? " with "+atoms[0] : ""));
					}
					if(collisionTime < aia.collisionTime()) {
						if (Debug.ON && Debug.DEBUG_NOW && (Debug.LEVEL > 2 || Debug.anyAtom(atoms))) {
							System.out.println("setting up time "+collisionTime+" for atom "+atoms[1]+(atoms.length > 1 ? " with "+atoms[0] : ""));
						}
						aia.setCollision(collisionTime, atoms[0], pHard);
					}//end if
				}//end if
			}
		}//end of actionPerformed
	} //end of collisionHandlerDown

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
    public static class Agent {  //need public so to use with instanceof
    	static PotentialHard nullPotential;
        public Atom atom, collisionPartner;
        public double collisionTime = Double.MAX_VALUE; //time to next collision
        public PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        public boolean movedInList = false;
        
        public Agent(Atom a) {atom = a;}
        
        public String toString() {
            return "Collider: "+atom.toString()+"; Partner: "+collisionPartner.toString()+"; Potential: "+collisionPotential.toString();
        }
        public final Atom atom() {return atom;}
        public final Atom collisionPartner() {return collisionPartner;}
        
        public void resetCollision() {
            collisionTime = periodCollisionTime();
            collisionPotential = nullPotential;
            collisionPartner = null;
        }
        
        //time to "collision" to update colliders for periodic boundaries
        protected double periodCollisionTime() {
            Space.Boundary boundary = atom.node.parentPhase().boundary();
            if(boundary instanceof Space.Boundary.Periodic) {
                if(!(atom.type instanceof AtomType.Sphere)) {return Double.MAX_VALUE;}
                Space.Vector p = atom.coord.momentum();
                Space.Vector dim = boundary.dimensions();
                double tmin = Double.MAX_VALUE;
                double d2 = 2.0*((AtomType.Sphere)atom.type).diameter(atom);
                int D = dim.D();
                for(int i=0; i<D; i++) {
                    double t = (dim.x(i)-d2)/p.x(i);
                    t = (t < 0) ? -t : t;//abs
                    tmin = (t < tmin) ? t : tmin;
                }
                return 0.25*atom.coord.mass()*tmin; //0.5*m*min of (dim.x/p.x, dim.y/p.y, etc.)
          //      return 0.25*atom.mass()*dim.D(p).abs().min(); //0.5*m*min of (dim.x/p.x, dim.y/p.y, etc.)
                //assumes range of potential is .le. diameter, simulation box is square (or x is smaller dimension)
            }
            return Double.MAX_VALUE;
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

