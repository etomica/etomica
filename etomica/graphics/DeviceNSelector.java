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
public class DeviceNSelector extends DeviceSlider {
    
    public DeviceNSelector(Simulation simulation, SpeciesAgent agent) {
        this(simulation, agent, new SimulationRestart(simulation));
    }
    
    public DeviceNSelector(Simulation simulation, SpeciesAgent agent, Action resetAction) {
        super(simulation.getController());
        this.resetAction = resetAction;
        
//        setNMajor(6);
	    setMinimum(0);
        int max = 60;
        if (agent.getNMolecules() > max) max = agent.getNMolecules();
	    setMaximum(max);
	    getSlider().setSnapToTicks(false);
	    getSlider().setMajorTickSpacing(10);
	    graphic(null).setSize(new java.awt.Dimension(40,30));
        setModifier(new ModifierNMolecule(agent));
        targetAction = new ActionGroupSeries(new Action[]{modifyAction,resetAction});
	    
	    if(agent.type.getSpecies().getName() == "") {
	        setLabel("Number of molecules");
	    } else {
	        setLabel(agent.type.getSpecies().getName() + " molecules");
	    }
    }

    /**
     * Returns the action used to "reset" the simulation after changing the 
     * number of molecules, SimulationRestart by default.
     */
    public Action getResetAction() {
        return resetAction;
    }
    
    protected final Action resetAction;
    
    //main method to demonstrate and test class
    public static void main(String[] args) {
        
        HSMD2D sim = new HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.species.setName("Disk");

        new DeviceNSelector(sim,sim.phase.getAgent(sim.species));

        graphic.makeAndDisplayFrame();
    }
    

} //end of DeviceNSelector
  