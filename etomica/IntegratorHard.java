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

    public String getVersion() {return "IntegratorHard:01.07.08/"+IntegratorHardAbstract.VERSION;}
    
    private AtomIterator upAtomIterator;
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
    private static final class CollisionHandlerUp implements Potential1Calculation, Potential2Calculation {
        double minCollisionTime;
        IntegratorHardAbstract.Agent aia;
        Atom atom1;
        
        public CollisionHandlerUp setAtom(Atom a) {
            atom1 = a;
            aia = (IntegratorHardAbstract.Agent)a.ia;
            minCollisionTime = aia.collisionTime();
            return this;
        }//end of setAtom
        
        //atom pair
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            Potential2Hard potentialHard = (Potential2Hard)potential;
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                if(pair.atom1() != atom1) setAtom(pair.atom1()); //need this if doing minimum collision time calculation for more than one atom
                double collisionTime = potentialHard.collisionTime(pair);
      /*debug          System.out.println("      UP "+pair.atom1.debugIndex+","
                                        +pair.atom2.debugIndex+","
                                        +(float)collisionTime+","
                                        +(float)minCollisionTime);*/
                if(collisionTime < minCollisionTime) {
                    minCollisionTime = collisionTime;
                    aia.setCollision(collisionTime, pair.atom2(), potentialHard);
                }//end if
            }//end while
        }//end of calculate(AtomPair...
        
        //single atom
        public void calculate(AtomIterator iterator, Potential1 potential) {
            if(!(potential instanceof Potential1Hard)) return;
            Potential1Hard potentialHard = (Potential1Hard)potential;
            while(iterator.hasNext()) {
                Atom atom = iterator.next();
                setAtom(atom);
                double collisionTime = potentialHard.collisionTime(atom);
                if(collisionTime < minCollisionTime) {
                    minCollisionTime = collisionTime;
                    aia.setCollision(collisionTime, null, potentialHard);
                }
            }//end while
        }//end of calculate(Atom...
    } //end of collisionHandlerUp

    //the down-handler has the logic of the Allen & Tildesley downList subroutine
    //sets collision times of atoms downlist of given atom to minimum of their current
    //value and their value with given atom
    private static final class CollisionHandlerDown implements Potential2Calculation {
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            Potential2Hard potentialHard = (Potential2Hard)potential;
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                double collisionTime = potentialHard.collisionTime(pair);
     /*debug           System.out.println("      DN "+pair.atom1.debugIndex+","
                                        +pair.atom2.debugIndex+","
                                        +(float)collisionTime+","
                                        +(float)((IntegratorHardAbstract.Agent)pair.atom2().ia).collisionTime());*/
                if(collisionTime < Double.MAX_VALUE) {
                    IntegratorHardAbstract.Agent aia = (IntegratorHardAbstract.Agent)pair.atom2().ia;
                    if(collisionTime < aia.collisionTime()) {
                        aia.setCollision(collisionTime, pair.atom1(), potentialHard);
                    }//end if
                }//end if
            }//end while
        }//end of calculate
    } //end of collisionHandlerDown

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
            if(aPartner != null && (aPartner == partner || aPartner == collider)) {//aPartner != null handles case where aPartner and partner are both null
                aAgent.resetCollision();
                potential.calculate(upList.set(a), collisionHandlerUp.setAtom(a));
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
                    potential.findCollisions(atom, UP, collisionHandlerUp);
                }
            }*/
            //to keep collision lists perfect, should do an upList on atoms that had this
            //atom on its neighbor list, but no longer do because it has moved away


        colliderAgent.resetCollision();
        potential.calculate(upList.set(collider), collisionHandlerUp.setAtom(collider));
        potential.calculate(downList.set(collider), collisionHandlerDown);
        if(partner != null) {
            ((IntegratorHardAbstract.Agent)partner.ia).resetCollision();
            potential.calculate(upList.set(partner), collisionHandlerUp.setAtom(partner));
            potential.calculate(downList.set(partner), collisionHandlerDown);
        }

    }//end of processCollision
    
    protected void updateAtom(Atom a) {
        Agent agent = (Agent)a.ia;
        agent.resetCollision();
        potential.calculate(upList.set(a), collisionHandlerUp.setAtom(a));
        potential.calculate(downList.set(a), collisionHandlerDown);
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
        if(isothermal) scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            ((IntegratorHardAbstract.Agent)atomIterator.next().ia).resetCollision();
        }
        potential.set(firstPhase).calculate(upList.set(), collisionHandlerUp); //assumes only one phase
        findNextCollider();
    }
                
    /**
    * Produces the Agent defined by this integrator.
    * One instance of an Agent is placed in each atom controlled by this integrator.
    */
        public Integrator.Agent makeAgent(Atom a) {
            return new IntegratorHardAbstract.Agent(a);
        }
             
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Simulation.instance = new Simulation(new Space2D());
	    IntegratorHard integratorHard1 = new IntegratorHard();
	 //   integratorHard1.setTimeStep(0.02);
	    SpeciesSpheres speciesSpheres1 = new SpeciesSpheres(10);
	    SpeciesSpheres speciesSpheres2 = new SpeciesSpheres(1);
	    speciesSpheres2.setColor(java.awt.Color.red);
	    final Phase phase = new Phase();
	    P2HardSphere potential12 = new P2HardSphere(4.0);
	    P2HardSphere potential22 = new P2HardSphere(5.0);
	    P2HardSphere potential11 = new P2HardSphere(3.0);
	    speciesSpheres2.setDiameter(5.0);
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
	    MeterRDF meterRDF = new MeterRDF();
	    DisplayPlot plot = new DisplayPlot();
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
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
                */
	//    displayPhase1.setColorScheme(integratorHard1.new HighlightColliders());
	    
	    displayPhase1.addDrawable(new Drawable() {
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
	                g.drawString(Integer.toString(a.index()), xP-20, yP-20);
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
	    });*/
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main

}//end of IntegratorHard

