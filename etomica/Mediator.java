package etomica;
import etomica.utility.HashMap2;

//Java2 imports
//import java.util.HashMap;
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.HashMap;
import etomica.utility.Iterator;
    
/** Class to perform actions that tie together the elements of a simulation.
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
        addMediatorPair(new Potential2Species.Default(this));
        addMediatorPair(new Potential1Species.Default(this));
        addMediatorPair(new PhaseSpecies.Default(this));
        addMediatorPair(new IntegratorPhase.Default(this));
        addMediatorPair(new DisplayIntegrator.Default(this));
        addMediatorPair(new Mediator.DisplayPhase.Default(this));
        addMediatorPair(new MeterPhase.Default(this));
        addMediatorPair(new MeterSpecies.Default(this));
        addMediatorPair(new ControllerIntegrator.Default(this));
        addMediatorPair(new DeviceNull.Default(this));
        addMediatorPair(new DisplayNull.Default(this));
        addMediatorPair(new DisplayMeter.Default(this));
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
            add((Simulation.Element)iter.next());
        }
    }
    
    /**
     * Processes the given element, tying it in with other appropriate elements that were added previously.
     * Loops through all submediators that are set to take a class of the given element (as identified by
     * its <code>baseClass</code> method, and allows each to do its processing of it.
     * Upon completion, the <code>setAdded</code> flag is set to <code>true</code>, so if the element
     * is added again, no action is taken.
     */
    public void add(Simulation.Element element) {
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
        protected Mediator mediator;
        public Subset(Mediator m) {mediator = m;}
        public abstract void add(Simulation.Element element);
        public abstract Class[] elementClasses();
        public void setPrior(Subset prior) {priorSubset = prior;}
        public final void setSuperceding(boolean b) {superceding = b;}
        public final boolean isSuperceding() {return superceding;}
    }//end of Subset
    
    public abstract static class Potential2Species extends Subset {
        public Potential2Species(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Potential2.class, Species.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Potential2) add((Potential2)element);
            if(element instanceof Species) add((Species)element);
        }
        
        public abstract void add(Potential2 p2);
        public abstract void add(Species s);
        
        public static class Default extends Potential2Species {
            public Default(Mediator m) {super(m);}
            /**
            * Adds given potential to parent simulations potential2 array.
            * Uses potential's species index values to assign element.
            */
            public void add(Potential2 p2) {
                //assumes speciesIndex of p2 is in bounds of array
                //should move to approach in Potential2 that sets species index by passing species
                //instead of an integer
                Potential2[][] potential2 = mediator.parentSimulation().potential2;
                potential2[p2.getSpecies1Index()][p2.getSpecies2Index()] = p2;
                potential2[p2.getSpecies2Index()][p2.getSpecies1Index()] = p2;
            }
            /**
            * Performs no action.
            */
            public void add(Species s) {}
        }//end of Default
    }//end of Potential2Species
    
    public abstract static class Potential1Species extends Subset {
        public Potential1Species(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Potential1.class, Species.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Potential1) add((Potential1)element);
            if(element instanceof Species) add((Species)element);
        }
        
        public abstract void add(Potential1 p1);
        public abstract void add(Species s);
        
        public static class Default extends Potential1Species {
            public Default(Mediator m) {super(m);}
            /**
            * Adds given potential to parent simulation's potential1 array.
            * Uses potential's species index values to assign element.
            */
            public void add(Potential1 p1) {
                //assumes speciesIndex of p1 is in bounds of array
                //should move to approach in Potential1 that sets species index by passing species
                //instead of an integer
                Potential1[] potential1 = mediator.parentSimulation().potential1;
                potential1[p1.getSpeciesIndex()] = p1;
            }
            /**
            * Performs no action.
            */
            public void add(Species s) {}
        }//end of Default
    }//end of Potential1Species
    
    public abstract static class PhaseSpecies extends Subset {
        public PhaseSpecies(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Phase.class, Species.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Phase) add((Phase)element);
            if(element instanceof Species) add((Species)element);
        }
        
        public abstract void add(Phase p2);
        public abstract void add(Species s);
        
        public static class Default extends PhaseSpecies {
            public Default(Mediator m) {super(m);}
            public void add(Phase phase) {
                for(Iterator is=mediator.parentSimulation().speciesList.iterator(); is.hasNext(); ) {
                    Species species = (Species)is.next();
                    if(species.wasAdded()) phase.addSpecies(species.makeAgent(phase));
                }
            }
            public void add(Species species) {
                for(Iterator ip=mediator.parentSimulation().phaseList.iterator(); ip.hasNext(); ) {
                    Phase phase = (Phase)ip.next();
                    if(phase.wasAdded()) phase.addSpecies(species.makeAgent(phase));
                }
            }
        }//end of Default
    }//end of PhaseSpecies
    
    public abstract static class IntegratorPhase extends Subset {
        public IntegratorPhase(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Phase.class, Integrator.class};}
        
        public void add(Simulation.Element element) {
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
                for(Iterator is=mediator.parentSimulation().integratorList.iterator(); is.hasNext(); ) {
                    Integrator integrator = (Integrator)is.next();
                    if(integrator.wasAdded() && integrator.wantsPhase()) phase.setIntegrator(integrator);
                    break;
                }
            }
            /**
             * Sets this as the integrator for all phases that have no integrator.
             */
            public void add(Integrator integrator) {
                for(Iterator ip=mediator.parentSimulation().phaseList.iterator(); ip.hasNext(); ) {
                    Phase phase = (Phase)ip.next();
                    if(phase.wasAdded() && phase.integrator() == null) phase.setIntegrator(integrator);
                }
            }
        }//end of Default
    }//end of IntegratorPhase

    public abstract static class DeviceIntegrator extends Subset {
        public DeviceIntegrator(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Device.class, Integrator.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Device) add((Device)element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Device d);
        public abstract void add(Integrator i);
        
        //no default
    }//end of DeviceIntegrator
    
    public abstract static class DisplayIntegrator extends Subset {
        public DisplayIntegrator(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Display.class, Integrator.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Display) add((Display)element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Display d);
        public abstract void add(Integrator i);
        
        public static class Default extends DisplayIntegrator {
            private boolean firstIntegrator = true;
            public Default(Mediator m) {super(m);}
            /**
             * Sets display as interval listener of first integrator, if it exists
             */
            public void add(Display display) {
                for(Iterator ip=mediator.parentSimulation().integratorList.iterator(); ip.hasNext(); ) {
                    Integrator integrator = (Integrator)ip.next();
                    if(integrator.wasAdded())  {
                        integrator.addIntervalListener(display);
                        break;
                    }
                }
            }
            /**
             * Registers all existing Displays as interval listeners of the added integrator, if it is the first one being added.
             */
            public void add(Integrator integrator) {
                if(!firstIntegrator) return;
                firstIntegrator = false;
                for(Iterator ip=mediator.parentSimulation().displayList.iterator(); ip.hasNext(); ) {
                    Display display = (Display)ip.next();
                    if(display.wasAdded()) integrator.addIntervalListener(display);
                }
            }
        }//end of Default
    }//end of DisplayIntegrator
    public abstract static class DisplayPhase extends Subset {
        public DisplayPhase(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Display.class, Phase.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Display) add((Display)element);
            if(element instanceof Phase) add((Phase)element);
        }
        
        public abstract void add(Display d);
        public abstract void add(Phase i);
        
        public static class Default extends Mediator.DisplayPhase {
            private Phase lastPhaseAdded = null;
            public Default(Mediator m) {super(m);}
            /**
             * Sets display's phase to be the last phase added, if it exists.
             * If a display is a DisplayPhase, assigns it to first phase lacking one.
             */
            public void add(Display display) {
                if(!vetoConnection(display,lastPhaseAdded)) display.setPhase(lastPhaseAdded);
                else if(display instanceof etomica.DisplayPhase) {
                    for(Iterator ip=mediator.parentSimulation().phaseList.iterator(); ip.hasNext(); ) {
                        Phase phase = (Phase)ip.next();
                        if(!vetoConnection(display, phase)) {
                            display.setPhase(phase);
                            break;
                        }//end if
                    }//end of for
                }//end of else if
            }
            /**
             * Sets the given phase, if it is the first added, as the phase of all displays added so far.
             * Also positions it to be the phase for any subsequently added displays (until another Phase is added).
             */
            public void add(Phase phase) {
                if(lastPhaseAdded == null) { //none of the displays yet have a phase; make this one it for all
                    for(Iterator ip=mediator.parentSimulation().displayList.iterator(); ip.hasNext(); ) {
                        Display display = (Display)ip.next();
                        if(display.wasAdded() && !vetoConnection(display, phase)) display.setPhase(phase);
                    }
                }
                lastPhaseAdded = phase;
            }
            //returns true if the given phase is already assigned to an etomica.DisplayPhase, false otherwise
            private boolean vetoConnection(Display display, Phase phase) {
                if(display == null || phase == null) return true;
                //no problem if display is not an etomica.DisplayPhase
                if(!(display instanceof etomica.DisplayPhase)) return false;
                //otherwise make sure phase doesn't already have a DisplayPhase
                for(Iterator ip=mediator.parentSimulation().displayList.iterator(); ip.hasNext(); ) {
                    Display d = (Display)ip.next();
                    if(d instanceof etomica.DisplayPhase && d.getPhase() == phase) return true;
                }
                //allow connection
                return false;
            }
                
        }//end of Default
    }//end of DisplayPhase
    public abstract static class MeterPhase extends Subset {
        public MeterPhase(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {MeterAbstract.class, Phase.class};}
        
        public void add(Simulation.Element element) {
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
                for(Iterator ip=mediator.parentSimulation().meterList.iterator(); ip.hasNext(); ) {
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
        
        public void add(Simulation.Element element) {
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
        
        public void add(Simulation.Element element) {
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
                for(Iterator ip=mediator.parentSimulation().controllerList.iterator(); ip.hasNext(); ) {
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
                for(Iterator ip=mediator.parentSimulation().integratorList.iterator(); ip.hasNext(); ) {
                    Integrator integrator = (Integrator)ip.next();
                    if(integrator.wasAdded()) controller.add(integrator);
                }
            }
        }//end of Default
    }//end of ControllerIntegrator

//temporary means to make a single-object mediator (does things for a single simulation element, instead of pair)
    public abstract static class DisplayNull extends Subset {
        public DisplayNull(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Display.class, null};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Display) add((Display)element);
        }
        
        public abstract void add(Display d);
        
        public static class Default extends DisplayNull {
            private final java.awt.GridBagConstraints gbcBox = new java.awt.GridBagConstraints();
            public Default(Mediator m) {
                super(m);
                gbcBox.gridx = 0;
            }
            /**
             * Adds displays graphic to the simulation display pane
             */
            public void add(Display display) {
                final java.awt.Component component = display.graphic(null);
                if(component == null) return; //display is not graphic
                if(display instanceof DisplayBox) {
                    mediator.parentSimulation().displayBoxPanel.add(component, gbcBox);
                }
                else {
                    mediator.parentSimulation().displayPanel.add(display.getLabel(),component);
                    //add a listener to update the tab label if the name of the display changes
                    display.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
                        public void propertyChange(java.beans.PropertyChangeEvent evt) {
                            if(evt.getPropertyName().equals("label")) {
                                int idx = mediator.parentSimulation.displayPanel.indexOfComponent(component);
                                mediator.parentSimulation().displayPanel.setTitleAt(idx,evt.getNewValue().toString());
                            }
                        }
                    });
                }
            }
        }//end of Default
        
        /**
         * Adding an instance of this pair mediator will cause the simulation to not
         * act on the addition of display objects.  This is useful if it is desired to
         * place the display elements in the simulation frame "by hand", instead of using
         * the simple default behavior.
         */
        public static class NoAction extends DisplayNull {
            public NoAction(Mediator m) {
                super(m);
                setSuperceding(true);//causes all previously added mediators to be ignored
            }
            
            public void add(Display display) {}
        }//end of NoAction
    }//end of DisplayNull
    public abstract static class DeviceNull extends Subset {
        public DeviceNull(Mediator m) {
            super(m);
        }

        public Class[] elementClasses() {return new Class[] {Device.class, null};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Device) add((Device)element);
        }
        
        public abstract void add(Device d);
        
        public static class Default extends DeviceNull {
            private final java.awt.GridBagConstraints gbc = new java.awt.GridBagConstraints();
            public Default(Mediator m) {
                super(m);
                gbc.gridx = 0;
            }
            /**
             * Adds displays graphic to the simulation display pane
             */
            public void add(Device device) {
                java.awt.Component component = device.graphic(null);
                if(device instanceof DeviceTable) {
                    mediator.parentSimulation().displayPanel.add(component);
                }
                else {
                    mediator.parentSimulation().devicePanel.add(component,gbc);
                }
            }
        }//end of Default

        /**
         * Adding an instance of this pair mediator will cause the simulation to not
         * act on the addition of device objects.  This is useful if it is desired to
         * place the device elements in the simulation frame "by hand", instead of using
         * the simple default behavior.
         */
        public static class NoAction extends DeviceNull {
            public NoAction(Mediator m) {
                super(m);
                setSuperceding(true);//causes all previously added mediators to be ignored
            }
            
            public void add(Device device) {}
        }//end of NoAction
    }//end of DeviceNull
    
    public abstract static class IntegratorNull extends Subset {
        public IntegratorNull(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Integrator.class, null};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Integrator d);
        
        //no Default defined
    }//end of IntegratorNull

    public abstract static class DisplayMeter extends Subset {
        public DisplayMeter(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {MeterAbstract.class, Display.class};}
        
        public void add(Simulation.Element element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof MeterAbstract) add((MeterAbstract)element);
            if(element instanceof Display) add((Display)element);
        }
        
        public abstract void add(MeterAbstract m);
        public abstract void add(Display d);
        
        public static class Default extends DisplayMeter {
            public Default(Mediator m) {super(m);}
            /**
             * Adds meter to all displays needing or accepting a meter
             */
            public void add(MeterAbstract meter) {
                for(Iterator ip=mediator.parentSimulation().displayList.iterator(); ip.hasNext(); ) {
                    Display display = (Display)ip.next();
                    if(display.wasAdded()) connect(display, meter);
                }//end of for
            }//end of add(meter)
            /**
             * Sets the given phase, if it is the first added, as the phase of all meters.
             */
            public void add(Display display) {
                for(Iterator ip=mediator.parentSimulation().meterList.iterator(); ip.hasNext(); ) {
                    MeterAbstract meter = (MeterAbstract)ip.next();
                    if(meter.wasAdded()) connect(display, meter);
                    
                }
            }//end of add(display)
            
            //need to handle MeterGroup
            private void connect(Display display, MeterAbstract meter) {
                if(meter instanceof Meter) {
                    if(display instanceof DatumSource.MultiUser) {
                        ((DatumSource.MultiUser)display).addDatumSources((Meter)meter);
                    }
                    else if(display instanceof DatumSource.User && ((DatumSource.User)display).getDatumSource() == null) {
                        ((DatumSource.User)display).setDatumSource((Meter)meter);
                    }
                }
                else if(meter instanceof MeterFunction) {
                    if(display instanceof DataSource.MultiUser) {
                        ((DataSource.MultiUser)display).addDataSources((MeterFunction)meter);
                    }
                    else if(display instanceof DataSource.User && ((DataSource.User)display).getDataSource() == null) {
                        ((DataSource.User)display).setDataSource((MeterFunction)meter);
                    }
                }
            }//end of connect
        }//end of Default (DisplayMeter)
        /**
         * Adding an instance of this pair mediator will cause the simulation to not
         * act to connect display and meter objects.  This is useful if it is desired to
         * connect these elements "by hand", instead of using
         * the simple default behavior.
         */
        public static class NoAction extends DisplayMeter {
            public NoAction(Mediator m) {
                super(m);
                setSuperceding(true);//causes all previously added mediators to be ignored
            }
            
            public void add(Display display) {}
            public void add(MeterAbstract meter) {}
        }//end of NoAction (DisplayMeter)
    }//end of DisplayMeter
    
}//end of Mediator