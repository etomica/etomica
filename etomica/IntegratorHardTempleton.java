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
public class IntegratorHardTempleton extends IntegratorHardAbstract implements EtomicaElement {

	public String getVersion() {return "IntegratorHard:01.07.08/"+IntegratorHardAbstract.VERSION;}
    
	private static final IteratorDirective upList = new IteratorDirective(IteratorDirective.UP);
	private static final IteratorDirective downList = new IteratorDirective(IteratorDirective.DOWN);
	private final CollisionHandlerUp collisionHandlerUp = new CollisionHandlerUp();
	private final CollisionHandlerDown collisionHandlerDown = new CollisionHandlerDown();

	//collision handler is passed to the potential and is notified of each collision
	//the potential detects.  The collision handler contains the necessary logic to
	//process this information so that the collision lists are kept up to date.
	//The potential should call the handler's setPotential method, with itself as the
	//argument, before beginning to detect collisions. 
    
	//the up-handler has the logic of the Allen & Tildesley upList subroutine
	//sets collision time of given atom to minimum value for collisions with all atoms uplist of it
	private static final class CollisionHandlerUp extends PotentialCalculation {
		double minCollisionTime;
		IntegratorHardAbstract.Agent aia;
		Atom atom1;
		private Potential2.Hard p2Hard;
		private Potential1.Hard p1Hard;
        
		public CollisionHandlerUp setAtom(Atom a) {
			atom1 = a;
			aia = (IntegratorHardAbstract.Agent)a.ia;
			minCollisionTime = aia.collisionTime();
			return this;
		}//end of setAtom
        
		//atom pair
		public void actionPerformed(AtomPair pair) {
			if(pair.atom1() != atom1) setAtom(pair.atom1()); //need this if doing minimum collision time calculation for more than one atom
			double collisionTime = p2Hard.collisionTime(pair);
  /*debug* /   System.out.println("      UP "+pair.atom1.toString()+","
									+pair.atom2.toString()+","
									+(float)collisionTime+","
									+(float)minCollisionTime); /* */
			if(collisionTime < minCollisionTime) {
				minCollisionTime = collisionTime;
				aia.setCollision(collisionTime, pair.atom2(), p2Hard);
			}//end if
   //info         System.out.println(atom1.toString()+" UP: "+count);
		}//end of calculate(AtomPair...
        
		//single atom
		public void actionPerformed(Atom atom) {
			setAtom(atom);
			double collisionTime = p1Hard.collisionTime(atom);
			if(collisionTime < minCollisionTime) {
				minCollisionTime = collisionTime;
				aia.setCollision(collisionTime, null, p1Hard);
			}
		}//end of calculate(Atom...
		public PotentialCalculation set(Potential1 p1) {
			if(!(p1 instanceof Potential1.Hard)) throw new RuntimeException("Error: Only hard potentials can be used with IntegratorHard");
			p1Hard = (Potential1.Hard)p1;
			return super.set(p1);
		}
		public PotentialCalculation set(Potential2 p2) {
			if(!(p2 instanceof Potential2.Hard)) throw new RuntimeException("Error: Only hard potentials can be used with IntegratorHard");
			p2Hard = (Potential2.Hard)p2;
			return super.set(p2);
		}
	} //end of collisionHandlerUp

	//the down-handler has the logic of the Allen & Tildesley downList subroutine
	//sets collision times of atoms downlist of given atom to minimum of their current
	//value and their value with given atom
	private static final class CollisionHandlerDown extends PotentialCalculation {
		private Potential2.Hard p2Hard;
		private Potential1.Hard p1Hard;
		public void actionPerformed(AtomPair pair) {
			double collisionTime = p2Hard.collisionTime(pair);
 /*debug * /   System.out.println("      DN "+pair.atom1.toString()+","
									+pair.atom2.toString()+","
									+(float)collisionTime+","
									+(float)((IntegratorHardAbstract.Agent)pair.atom2().ia).collisionTime());  /* */
			if(collisionTime < Double.MAX_VALUE) {
				IntegratorHardAbstract.Agent aia = (IntegratorHardAbstract.Agent)pair.atom2().ia;
				if(collisionTime < aia.collisionTime()) {
					aia.setCollision(collisionTime, pair.atom1(), p2Hard);
				}//end if
			}//end if
//info            System.out.println(" DN: "+count);
		}//end of actionPerformed
		public PotentialCalculation set(Potential1 p1) {
			if(!(p1 instanceof Potential1.Hard)) throw new RuntimeException("Error: Only hard potentials can be used with IntegratorHard");
			p1Hard = (Potential1.Hard)p1;
			return super.set(p1);
		}
		public PotentialCalculation set(Potential2 p2) {
			if(!(p2 instanceof Potential2.Hard)) throw new RuntimeException("Error: Only hard potentials can be used with IntegratorHard");
			p2Hard = (Potential2.Hard)p2;
			return super.set(p2);
		}
	} //end of collisionHandlerDown

	public IntegratorHardTempleton() {
		this(Simulation.instance);
	}
	public IntegratorHardTempleton(Simulation sim) {
		super(sim);
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("Collision-based molecular dynamics simulation of hard potentials");
		return info;
	}
    
   /**
	* Loops through all atoms to identify the one with the smallest value of collisionTime
	* Collision time is obtained from the value stored in the Integrator.Agent from each atom.
	*/
	protected void findNextCollider() {
		//find next collision pair by looking for minimum collisionTime
		double minCollisionTime = Double.MAX_VALUE;
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			IntegratorHardAbstract.Agent ia = (IntegratorHardAbstract.Agent)atomIterator.next().ia;
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
			Atom a = atomIterator.next();
			if(a == collider) { //finished with atoms before collider...
				if(partner == null) break;  //and there is no partner, so we're done, or...
				else continue;              //...else just skip this atom and continue with loop
			}
			if(a == partner) break; //finished with atoms before partner; we're done
			IntegratorHardAbstract.Agent aAgent = (IntegratorHardAbstract.Agent)a.ia;
			Atom aPartner = aAgent.collisionPartner();
			if(aPartner != null && (aPartner == partner || aPartner == collider)) {//aPartner != null handles case where aPartner and partner are both null
				aAgent.resetCollision();
				potential.calculate(firstPhase, upList.set(a), collisionHandlerUp.setAtom(a));
			}
		}//end while
			//reset collision partners of atoms that are now up from this atom but still list it as their
			//collision partner.  Assumes this atom was moved down list, but this won't always be the case
			//This bit could be made more efficient
            
			//if(a movedInList) {  add a means for bump method to declare it moved atom in the list
	   /*     atomIterator.reset(a);
			while(atomIterator.hasNext()) {
				Atom atom = atomIterator.next();
				if(((Agent)atom.ia).collisionPartner == a) {  //upList atom could have atom as collision partner if atom was just moved down list
					potential.findCollisions(atom, UP, collisionHandlerUp);
				}
			}*/
			//to keep collision lists perfect, should do an upList on atoms that had this
			//atom on its neighbor list, but no longer do because it has moved away


		colliderAgent.resetCollision();
		potential.calculate(firstPhase, upList.set(collider), collisionHandlerUp.setAtom(collider));
		potential.calculate(firstPhase, downList.set(collider), collisionHandlerDown);
		if(partner != null) {
			((IntegratorHardAbstract.Agent)partner.ia).resetCollision();
			potential.calculate(firstPhase, upList.set(partner), collisionHandlerUp.setAtom(partner));
			potential.calculate(firstPhase, downList.set(partner), collisionHandlerDown);
		}

	}//end of processCollision
    
	protected void updateAtom(Atom a) {
		Agent agent = (Agent)a.ia;
		agent.resetCollision();
		potential.calculate(firstPhase, upList.set(a), collisionHandlerUp.setAtom(a));
		potential.calculate(firstPhase, downList.set(a), collisionHandlerDown);
	}

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
			if(a.coord.isStationary()) {continue;}  //skip if atom is stationary
			a.coord.freeFlight(tStep);
	//     a.translateBy(tStep*a.rm(),a.momentum());
		}
	}

	/**
	* Sets up the integrator for action.
	* Do an upList call for each atom and find the next collider
	*/
	protected void doReset() {
		if(isothermal) scaleMomenta(Math.sqrt(this.temperature/meterTemperature.currentValue(firstPhase.speciesMaster)));
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			((IntegratorHardAbstract.Agent)atomIterator.next().ia).resetCollision();
		}
		potential.calculate(firstPhase, upList.set(), collisionHandlerUp); //assumes only one phase
		findNextCollider();
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
					IntegratorHardAbstract.Agent agent = (IntegratorHardAbstract.Agent)a.ia;
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
			public void collisionAction(IntegratorHardAbstract.Agent agent) {
				int i0 = agent.atom.debugIndex;
				int i1 = (agent.collisionPartner!=null)?agent.collisionPartner.debugIndex:-1;
	  //          System.out.print(i0+" "+i1+" ");
				for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {
					IntegratorHardAbstract.Agent aAgent = (IntegratorHardAbstract.Agent)a.ia;
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

