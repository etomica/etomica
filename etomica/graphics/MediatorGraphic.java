package etomica.graphics;
import etomica.*;

//Java2 imports
//import java.util.HashMap;
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.java2.Iterator;
    
/** 
 * Class to perform actions that tie together the elements of a simulation.
 */
    
public class MediatorGraphic extends Mediator {

    public static String getVersion() {return "MediatorGraphic:01.11.20";}
    private SimulationGraphic parentSimulation;
    
    public MediatorGraphic(SimulationGraphic sim) {
        super(sim);
        parentSimulation = sim;
        addMediatorPair(new DisplayIntegrator.Default(this));
        addMediatorPair(new MediatorGraphic.DisplayPhase.Default(this));
        addMediatorPair(new DeviceNull.Default(this));
        addMediatorPair(new DisplayNull.Default(this));
        addMediatorPair(new DisplayMeter.Default(this));
        addMediatorPair(new ControllerNullDefault(this));
    }
        
    public abstract static class DeviceIntegrator extends Mediator.Subset {
        public DeviceIntegrator(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Device.class, Integrator.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Device) add((Device)element);
            if(element instanceof Integrator) add((Integrator)element);
        }
        
        public abstract void add(Device d);
        public abstract void add(Integrator i);
        
        //no default
    }//end of DeviceIntegrator
    
    public abstract static class DisplayIntegrator extends Mediator.Subset {
        public DisplayIntegrator(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Display.class, Integrator.class};}
        
        public void add(SimulationElement element) {
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
                for(Iterator ip=mediator.parentSimulation().integratorList().iterator(); ip.hasNext(); ) {
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
                for(Iterator ip=((SimulationGraphic)mediator.parentSimulation()).displayList().iterator(); ip.hasNext(); ) {
                    Display display = (Display)ip.next();
                    if(display.wasAdded()) integrator.addIntervalListener(display);
                }
            }
        }//end of Default
    }//end of DisplayIntegrator
    public abstract static class DisplayPhase extends Mediator.Subset {
        public DisplayPhase(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Display.class, Phase.class};}
        
        public void add(SimulationElement element) {
            if(!superceding && priorSubset != null) priorSubset.add(element);
            if(element instanceof Display) add((Display)element);
            if(element instanceof Phase) add((Phase)element);
        }
        
        public abstract void add(Display d);
        public abstract void add(Phase i);
        
        public static class Default extends MediatorGraphic.DisplayPhase {
            private Phase lastPhaseAdded = null;
            public Default(Mediator m) {super(m);}
            /**
             * Sets display's phase to be the last phase added, if it exists.
             * If a display is a DisplayPhase, assigns it to first phase lacking one.
             */
            public void add(Display display) {
                if(!vetoConnection(display,lastPhaseAdded)) display.setPhase(lastPhaseAdded);
                else if(display instanceof etomica.graphics.DisplayPhase) {
                    for(Iterator ip=mediator.parentSimulation().phaseList().iterator(); ip.hasNext(); ) {
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
                    for(Iterator ip=((SimulationGraphic)mediator.parentSimulation()).displayList().iterator(); ip.hasNext(); ) {
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
                if(!(display instanceof etomica.graphics.DisplayPhase)) return false;
                //otherwise make sure phase doesn't already have a DisplayPhase
                for(Iterator ip=((SimulationGraphic)mediator.parentSimulation()).displayList().iterator(); ip.hasNext(); ) {
                    Display d = (Display)ip.next();
                    if(d instanceof etomica.graphics.DisplayPhase && d.getPhase() == phase) return true;
                }
                //allow connection
                return false;
            }
                
        }//end of Default
    }//end of DisplayPhase


//temporary means to make a single-object mediator (does things for a single simulation element, instead of pair)
    public abstract static class DisplayNull extends Mediator.Subset {
        public DisplayNull(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {Display.class, null};}
        
        public void add(SimulationElement element) {
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
                    ((SimulationGraphic)mediator.parentSimulation()).panel().displayBoxPanel.add(component, gbcBox);
                }
                else {
                    ((SimulationGraphic)mediator.parentSimulation()).panel().displayPanel.add(display.getLabel(),component);
                    //add a listener to update the tab label if the name of the display changes
                    display.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
                        public void propertyChange(java.beans.PropertyChangeEvent evt) {
                            if(evt.getPropertyName().equals("label")) {
                                int idx = ((SimulationGraphic)mediator.parentSimulation()).panel().displayPanel.indexOfComponent(component);
                                ((SimulationGraphic)mediator.parentSimulation()).panel().displayPanel.setTitleAt(idx,evt.getNewValue().toString());
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
    public abstract static class DeviceNull extends Mediator.Subset {
        public DeviceNull(Mediator m) {
            super(m);
        }

        public Class[] elementClasses() {return new Class[] {Device.class, null};}
        
        public void add(SimulationElement element) {
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
                    ((SimulationGraphic)mediator.parentSimulation()).panel().displayPanel.add(component);
                }
                else {
                    ((SimulationGraphic)mediator.parentSimulation()).panel().devicePanel.add(component,gbc);
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
    
    public abstract static class DisplayMeter extends Mediator.Subset {
        public DisplayMeter(Mediator m) {super(m);}

        public Class[] elementClasses() {return new Class[] {MeterAbstract.class, Display.class};}
        
        public void add(SimulationElement element) {
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
                for(Iterator ip=((SimulationGraphic)mediator.parentSimulation()).displayList().iterator(); ip.hasNext(); ) {
                    Display display = (Display)ip.next();
                    if(display.wasAdded()) connect(display, meter);
                }//end of for
            }//end of add(meter)
            /**
             * Sets the given phase, if it is the first added, as the phase of all meters.
             */
            public void add(Display display) {
                for(Iterator ip=mediator.parentSimulation().meterList().iterator(); ip.hasNext(); ) {
                    MeterAbstract meter = (MeterAbstract)ip.next();
                    if(meter.wasAdded()) connect(display, meter);
                    
                }
            }//end of add(display)
            
            //need to handle MeterGroup
            private void connect(Display display, MeterAbstract meter) {
                if(meter instanceof MeterScalar) {
                    if(display instanceof DatumSource.MultiUser) {
                        ((DatumSource.MultiUser)display).addDatumSources((MeterScalar)meter);
                    }
                    else if(display instanceof DatumSource.User && ((DatumSource.User)display).getDatumSource() == null) {
                        ((DatumSource.User)display).setDatumSource((MeterScalar)meter);
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
    
    /**
     * Causes addition of button that toggles controller state.
     */
    public static class ControllerNullDefault extends ControllerNull {
        public ControllerNullDefault(Mediator m) {
            super(m);
        }
        /**
         * Causes addition of button that toggles controller state.
         */
        public void add(Controller controller) {
            mediator.add(new DeviceControllerButton(mediator.parentSimulation(), controller));
            
        }
    }//end of ControllerNullDefault

    
}//end of MediatorGraphic