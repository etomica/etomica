package etomica;
import etomica.utility.HashMap2;

//Java2 imports
//import java.util.HashMap;
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.java2.HashMap;
import etomica.utility.java2.Iterator;
    
/** 
 * Class to perform actions that tie together the elements of a simulation.
 * Such actions include placement of species agents in phases, and registering meters and displays with
 * the integrator.  
 * Implemented by calling the go() method of the class before starting the simulation. Only the first call of
 * this method performs any action; subsequent calls are ignored, so the method should be invoked only
 * after all components are in place.<br>
 * This facility is provided mainly as a convenience to visual programming environments.
 * All of the actions performed by the ElementCoordinator may instead be done directly in the applet
 * or application.  If multiple phases or integrators are present in a simulation, it may be necessary
 * to coordinate them directly rather than invoking an element coordinator.<br>
 * A default element coordinator is held in the Simulation class.  It is invoked automatically by the
 * start() and run() methods of the Controller class.  If a DisplayPhase object is being used, the
 * element coordinator should be invoked directly by the applet or application to ensure that the initial
 * configuration is displayed (i.e., so that the configuration shows up before the controller is started).
 * To override the choice of the default element coordinator, be sure to change the handle assigned
 * to the Default.elementCoordinator field.  To have the ElementCoordinator take no action
 * when its go method is called, set its <code>completed</code> field to <code>true</code>.
 */
    
public class Mediator implements java.io.Serializable {

    public static String getVersion() {return "01.03.11";}
    private Simulation parentSimulation;
    private final HashMap2 mediatorTable = new HashMap2();
    
    public Mediator(Simulation sim) {
        parentSimulation = sim;
        addMediatorPair(new PhaseSpecies.Default(this));
        addMediatorPair(new IntegratorPhase.Default(this));
        addMediatorPair(new MeterPhase.Default(this));
        addMediatorPair(new MeterSpecies.Default(this));
        addMediatorPair(new ControllerIntegrator.Default(this));
        addMediatorPair(new PhasePotential.Default(this));
        addMediatorPair(new PotentialSpecies.Default(this));
    }
    public Simulation parentSimulation() {return parentSimulation;}
    public Class baseClass() {return Mediator.class;}
    
    public void addMediatorPair(Subset mPair) {
        Class[] classes = mPair.elementClasses();
        if(classes.length != 2) return;  //should throw exception
        //registers a pair mediator in the pair hashtable.  If another entry was there,
        //it is returned and held in the "prior" field of the new one.
        mPair.setPrior((Subset)mediatorTable.put(classes[0], classes[1], mPair));
    }
        
    /**
    * Flag indicating whether the element coordinator has yet been invoked.
    */
    public boolean completed = false;
    /**
    * Processes all simulation elements through the mediator.  
    * This method should be called only if elements were not added individually
    * via calls to mediators add method.  Method has effect only first time
    * it is called.  Subsequent calls do nothing, unless <code>completed</code>
    * field is set to <code>false</code>.
    */
    public void go() {
 //       if(completed) return;
        completed = true;
        for(Iterator iter=parentSimulation().allElements().iterator(); iter.hasNext();) {
            add((SimulationElement)iter.next());
        }
    }
    
    /**
     * Processes the given element, tying it in with other appropriate elements that were added previously.
     * Loops through all submediators that are set to take a class of the given element (as identified by
     * its <code>baseClass</code> method, and allows each to do its processing of it.
     * Upon completion, the <code>setAdded</code> flag is set to <code>true</code>, so if the element
     * is added again, no action is taken.
     */
    public void add(SimulationElement element) {
        if(element.wasAdded()) return;
        HashMap table = mediatorTable.get(element.baseClass());
        if(table == null) return;
        for(Iterator iter=table.values().iterator(); iter.hasNext(); ) {
            Subset subMediator = (Subset)iter.next();
            subMediator.add(element);
        }
        element.setAdded(true);
    }
    
    /**
     * A Mediator specifically designed to process a subset of all possible classes (usually two of them).
     * The classes processed are returned by the <code>elementClasses</code> method.
     */
    public abstract static class Subset implements java.io.Serializable {
        protected Subset priorSubset;
        protected boolean superceding = false;
        public final Mediator mediator;
        public Subset(Mediator m) {mediator = m;}
        public abstract void add(SimulationElement element);
        public abstract Class[] elementClasses();
        public void setPrior(Subset prior) {priorSubset = prior;}
        public final void setSuperceding(boolean b) {superceding = b;}
        public final boolean isSuperceding() {return superceding;}
    }//end of Subset
    
    
    public abstract static class PhaseSpecies extends Subset {
        public PhaseSpecies(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Phase.class, Species.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Phase) add((Phase)element);
            if(element instanceof Species) add((Species)element);
        }
        
        public abstract void add(Phase p2);
        public abstract void add(Species s);
        
        public static class Default extends PhaseSpecies {
            public Default(Mediator m) {super(m);}
            public void add(Phase phase) {
                for(Iterator is=mediator.parentSimulation().speciesList().iterator(); is.hasNext(); ) {
                    Species species = (Species)is.next();
                    if(species.wasAdded()) phase.speciesMaster.addSpecies(species);
                }
            }
            public void add(Species species) {
                for(Iterator ip=mediator.parentSimulation().phaseList().iterator(); ip.hasNext(); ) {
                    Phase phase = (Phase)ip.next();
                    if(phase.wasAdded()) phase.speciesMaster.addSpecies(species);
                }
            }
        }//end of Default
    }//end of PhaseSpecies
    
   public abstract static class PhasePotential extends Subset {
        public PhasePotential(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Phase.class, Potential.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Phase) add((Phase)element);
            if(element instanceof Potential) add((Potential)element);
        }
        
        public abstract void add(Phase phase);
        public abstract void add(Potential potential);
        
        //does nothing
        public static class Default extends PhasePotential {
            public Default(Mediator m) {super(m);}
            public void add(Phase phase) {
                for(Iterator is=mediator.parentSimulation().potentialList().iterator(); is.hasNext(); ) {
                    Potential potential = (Potential)is.next();
//                    if(potential.wasAdded()) phase.addPotential(potential);
                }
            }
            public void add(Potential potential) {
                for(Iterator ip=mediator.parentSimulation().phaseList().iterator(); ip.hasNext(); ) {
                    Phase phase = (Phase)ip.next();
//                    if(phase.wasAdded()) phase.addPotential(potential);
                }
            }
        }//end of Default
    }//end of PhasePotential
    

    public abstract static class PhaseNull extends Subset {
        public PhaseNull(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Phase.class, null};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Phase) add((Phase)element);
        }
        
        public abstract void add(Phase d);
        
        //no Default defined
    }//end of PhaseNull

    public abstract static class PotentialSpecies extends Subset {
        public PotentialSpecies(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Potential.class, Species.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Species) add((Species)element);
            if(element instanceof Potential) add((Potential)element);
        }
        
        public abstract void add(Species species);
        public abstract void add(Potential potential);
        
        //very crude mediator -- assigns species to any potentials needing them.
        //probably works correctly only if there is one species
        public static class Default extends PotentialSpecies {
            public Default(Mediator m) {super(m);}
            //set for all potentials needing a species
            public void add(Species species) {
                for(Iterator is=mediator.parentSimulation().potentialList().iterator(); is.hasNext(); ) {
                    Potential potential = (Potential)is.next();
                    if(potential.wasAdded() 
                        && potential.getSpecies() == null 
                        && potential.parentPotential() instanceof PotentialMaster) 
                            potential.setSpecies(new Species[] {species});
                }
            }
            public void add(Potential potential) {
                //do nothing if potential already has species assigned, or if it is not a molecule potential
                if(!(potential.getSpecies() == null && potential.parentPotential() instanceof PotentialMaster)) return;
                //set using first species
                for(Iterator ip=mediator.parentSimulation().speciesList().iterator(); ip.hasNext(); ) {
                    Species species = (Species)ip.next();
                    if(species.wasAdded()) {
                        potential.setSpecies(new Species[] {species});
                        return;
                    }
                }
            }
        }//end of Default
    }//end of PhasePotential

       public abstract static class IntegratorPhase extends Subset {
        public IntegratorPhase(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Phase.class, Integrator.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Phase) add((Phase)element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Phase p);
        public abstract void add(Integrator i);
        
        public static class Default extends IntegratorPhase {
            public Default(Mediator m) {super(m);}
            /**
             * Loops over all integrators and adds phase to first that wants a phase
             */
            public void add(Phase phase) {
                for(Iterator is=mediator.parentSimulation().integratorList().iterator(); is.hasNext(); ) {
                    Integrator integrator = (Integrator)is.next();
                    if(integrator.wasAdded() && integrator.wantsPhase() && phase.integrator()==null) {
                        phase.setIntegrator(integrator);
                        integrator.addIntervalListener(new PhaseAction.ImposePbc(phase));
                    }
                    break;
                }
            }
            /**
             * Sets this as the integrator for all phases that have no integrator.
             */
            public void add(Integrator integrator) {
                for(Iterator ip=mediator.parentSimulation().phaseList().iterator(); ip.hasNext(); ) {
                    Phase phase = (Phase)ip.next();
                    if(phase.wasAdded() && integrator.wantsPhase() && phase.integrator() == null) { 
                        phase.setIntegrator(integrator);
                        integrator.addIntervalListener(new PhaseAction.ImposePbc(phase));
                    }
                }
            }
        }//end of Default
        
        public static class NoCentralImage extends IntegratorPhase {
            public NoCentralImage(Mediator m) {
                super(m);
                setSuperceding(true);
            }
            /**
             * Loops over all integrators and adds phase to first that wants a phase
             */
            public void add(Phase phase) {
                for(Iterator is=mediator.parentSimulation().integratorList().iterator(); is.hasNext(); ) {
                    Integrator integrator = (Integrator)is.next();
                    if(integrator.wasAdded() && integrator.wantsPhase() && phase.integrator()==null) {
                        phase.setIntegrator(integrator);
                    }
                    break;
                }
            }
            /**
             * Sets this as the integrator for all phases that have no integrator.
             */
            public void add(Integrator integrator) {
                for(Iterator ip=mediator.parentSimulation().phaseList().iterator(); ip.hasNext(); ) {
                    Phase phase = (Phase)ip.next();
                    if(phase.wasAdded() && integrator.wantsPhase() && phase.integrator() == null) { 
                        phase.setIntegrator(integrator);
                    }
                }
            }
        }//end of NoCentralImage
        
    }//end of IntegratorPhase

    public abstract static class MeterPhase extends Subset {
        public MeterPhase(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {MeterAbstract.class, Phase.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof MeterAbstract) add((MeterAbstract)element);
            if(element instanceof Phase) add((Phase)element);
        }
        
        public abstract void add(MeterAbstract m);
        public abstract void add(Phase p);
        
        public static class Default extends MeterPhase {
            private Phase lastPhaseAdded = null;
            public Default(Mediator m) {super(m);}
            /**
             * Sets meter's phase to be the last phase added, if it exists, and if 
             * no phase has been previously set for this meter.
             */
            public void add(MeterAbstract meter) {
                if(lastPhaseAdded != null && meter.getPhase() == null) meter.setPhase(lastPhaseAdded);
            }
            /**
             * Sets the given phase as the phase of all previously added meters, and positions
             * it to be the phase of all subsequently added meters, until another phase is added.
             */
            public void add(Phase phase) {
                for(Iterator ip=mediator.parentSimulation().meterList().iterator(); ip.hasNext(); ) {
                    MeterAbstract meter = (MeterAbstract)ip.next();
                    if(meter.wasAdded() && meter.getPhase() == null) meter.setPhase(phase);
                }
                lastPhaseAdded = phase;
            }
        }//end of Default
    }//end of MeterPhase
    public abstract static class MeterSpecies extends Subset {
        public MeterSpecies(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {MeterAbstract.class, Species.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof MeterAbstract) add((MeterAbstract)element);
            if(element instanceof Species) add((Species)element);
        }
        
        public abstract void add(MeterAbstract m);
        public abstract void add(Species s);
        
        public static class Default extends MeterSpecies {
            public Default(Mediator m) {super(m);}
            /**
             * Does nothing.
             */
            public void add(MeterAbstract meter) {
            }
            /**
             * Does nothing
             */
            public void add(Species species) {
            }
        }//end of Default
    }//end of MeterSpecies

    public abstract static class ControllerIntegrator extends Subset {
        public ControllerIntegrator(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Controller.class, Integrator.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Controller) add((Controller)element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Controller d);
        public abstract void add(Integrator i);
        
        public static class Default extends ControllerIntegrator {
            private boolean firstController = true;
            public Default(Mediator m) {super(m);}
            /**
             * Sets first controller as controller of given integrator
             */
            public void add(Integrator integrator) {
                for(Iterator ip=mediator.parentSimulation().controllerList().iterator(); ip.hasNext(); ) {
                    Controller controller = (Controller)ip.next();
                    if(controller.wasAdded()) {
                        controller.add(integrator);
                        break;
                    }
                }
            }
            /**
             * Sets given controller, if it is the first added, as controller of all integrators.
             */
            public void add(Controller controller) {
                if(!firstController) return;
                firstController = false;
                for(Iterator ip=mediator.parentSimulation().integratorList().iterator(); ip.hasNext(); ) {
                    Integrator integrator = (Integrator)ip.next();
                    if(integrator.wasAdded()) controller.add(integrator);
                }
            }
        }//end of Default
    }//end of ControllerIntegrator

    public abstract static class IntegratorNull extends Subset {
        public IntegratorNull(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Integrator.class, null};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Integrator d);
        
        //no Default defined
    }//end of IntegratorNull

    public abstract static class ControllerNull extends Mediator.Subset {
        public ControllerNull(Mediator m) {
            super(m);
        }

        public Class[] elementClasses() {return new Class[] {Controller.class, null};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Controller) add((Controller)element);
        }
        
        public abstract void add(Controller d);
        
        /**
         * Adding an instance of this pair mediator will cause the simulation to not
         * act on the addition of controller objects.  
         */
        public static class NoAction extends ControllerNull {
            public NoAction(Mediator m) {
                super(m);
                setSuperceding(true);//causes all previously added mediators to be ignored
            }
            
            public void add(Controller controller) {}
        }//end of NoAction
    }//end of ControllerNull

}//end of Mediator