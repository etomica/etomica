package etomica.graphics;

import etomica.Controller;
import etomica.Integrator;
import etomica.Modifier;
import etomica.Phase;
import etomica.Simulation;
import etomica.SpeciesAgent;
import etomica.SpeciesSpheresMono;
import etomica.integrator.IntegratorHard;
import etomica.modifier.ModifierNMolecule;
import etomica.potential.P2HardSphere;
import etomica.simulations.HSMD2D;
import etomica.units.Dimension;

/**
 * Slider that selects the number of atoms of a given species in a phase.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/29/02 (DAK) new
  * 01/03/03 (DAK) changed snapToTicks to false
  */
  
  //could improve by allowing to set phase/agent before elementcoordinator call
public class DeviceNSelector extends DeviceSlider {
    
    private SpeciesAgent speciesAgent;
    private Integrator integrator;
    private etomica.action.SimulationRestart restartAction;
    private DisplayPhase display;
    
    public DeviceNSelector(Controller controller, SpeciesAgent agent) {
        super();
        this.speciesAgent = agent;
        this.integrator = integrator;
        restartAction = new etomica.action.SimulationRestart(sim);
        
        setModifier(new ModifierNMolecule());
//        setNMajor(6);
	    setMinimum(0);
	    setMaximum(60);
	    getSlider().setSnapToTicks(false);
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
    

       
    //main method to demonstrate and test class
    public static void main(String[] args) {
        
        HSMD2D sim = new HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.species.setName("Disk");
        
//        sim.elementCoordinator.go();
        DeviceNSelector nSelector = new DeviceNSelector(sim.phase.getAgent(sim.species));
        nSelector.setDisplayPhase(display);
        graphic.makeAndDisplayFrame();
    }//end of main
    

} //end of DeviceNSelector
  