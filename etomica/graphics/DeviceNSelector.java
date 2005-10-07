package etomica.graphics;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.atom.SpeciesAgent;
import etomica.modifier.ModifierNMolecule;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.HSMD2D;

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
    
    public DeviceNSelector(Simulation simulation, SpeciesAgent agent) {
        super(simulation.getController());
        SimulationRestart restartAction = new SimulationRestart(simulation);
        
//        setNMajor(6);
	    setMinimum(0);
        int max = 60;
        if (agent.getNMolecules() > max) max = agent.getNMolecules();
	    setMaximum(max);
	    getSlider().setSnapToTicks(false);
	    getSlider().setMajorTickSpacing(10);
	    graphic(null).setSize(new java.awt.Dimension(40,30));
        setModifier(new ModifierNMolecule(agent));
        targetAction = new ActionGroupSeries(new Action[]{modifyAction,restartAction});
	    
	    if(agent.type.getSpecies().getName() == "") {
	        setLabel("Number of molecules");
	    } else {
	        setLabel(agent.type.getSpecies().getName() + " molecules");
	    }
//	    setShowBorder(true);
    }
       
    //main method to demonstrate and test class
    public static void main(String[] args) {
        
        HSMD2D sim = new HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.species.setName("Disk");
        
//        sim.elementCoordinator.go();
        new DeviceNSelector(sim,sim.phase.getAgent(sim.species));
//        nSelector.setDisplayPhase(graphic.display);
        graphic.makeAndDisplayFrame();
    }//end of main
    

} //end of DeviceNSelector
  