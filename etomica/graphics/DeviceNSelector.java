package etomica.graphics;

import etomica.*;

/**
 * Slider that selects the number of atoms of a given species in a phase.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/29/02 (DAK) new
  */
  
  //could improve by allowing to set phase/agent before elementcoordinator call
public class DeviceNSelector extends DeviceSlider {
    
    private SpeciesAgent speciesAgent;
    private Integrator integrator;
    private etomica.action.SimulationRestart restartAction;
    private DisplayPhase display;
    
    public DeviceNSelector(SpeciesAgent agent) {
        super(agent.node.parentSimulation());
        this.speciesAgent = agent;
        this.integrator = agent.node.parentPhase().integrator();
        restartAction = new etomica.action.SimulationRestart(agent.node.parentSimulation());
        
        setModulator(new NMoleculeModulator());
//        setNMajor(6);
	    setMinimum(0);
	    setMaximum(60);
	    getSlider().setSnapToTicks(true);
	    getSlider().setMajorTickSpacing(10);
	    graphic(null).setSize(new java.awt.Dimension(40,30));
	    
	    if(agent.node.parentSpecies().getName() == "") {
	        setLabel("Number of molecules");
	    } else {
	        setLabel(agent.node.parentSpecies().getName() + " molecules");
	    }
//	    setShowBorder(true);
    }
    
    public void setDisplayPhase(DisplayPhase display) {
        this.display = display;
    }
    public Display getDisplayPhase() {return display;}
    

       
    //Modulator class to connect slider to species agents to change the
    //number of molecules
    private class NMoleculeModulator extends ModulatorAbstract {
        
        public NMoleculeModulator() {
            super();
        }
        
        public void setValue(double d) {
 //           if(initializing) return;
            if(d < 0) d = 0;
            boolean isPaused = integrator.isPaused();
   //         if(!isPaused) {
                try {
                    restartAction.actionPerformed();
                } catch (NullPointerException ex) {return;}
   //         }
            try {
                 speciesAgent.setNMolecules((int)d);
            } catch(NullPointerException ex) {}
                try {
                    restartAction.actionPerformed();
                } catch (NullPointerException ex) {return;}
            if(display != null) display.repaint();
            integrator.reset();
        }
        public double getValue() {return (speciesAgent!=null)?(double)speciesAgent.moleculeCount():0;}
        public etomica.units.Dimension getDimension() {return etomica.units.Dimension.NULL;}
    }//end of NMoleculeModulator
    
    //main method to demonstrate and test class
    public static void main(String[] args) {
        
        SimulationGraphic sim = new SimulationGraphic(new Space2D());
        Simulation.instance = sim;
        SpeciesSpheresMono species = new SpeciesSpheresMono();
        species.setName("Disk");
        Phase phase = new Phase();
        P2HardSphere potential = new P2HardSphere();
        IntegratorHard integrator = new IntegratorHard();
        Controller controller = new Controller();
        DisplayPhase display = new DisplayPhase();
        
        sim.elementCoordinator.go();
        DeviceNSelector nSelector = new DeviceNSelector(phase.getAgent(species));
        nSelector.setDisplayPhase(display);
        sim.elementCoordinator.go();
        sim.makeAndDisplayFrame();
    }//end of main
    

} //end of DeviceNSelector
  