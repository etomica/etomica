package etomica.graphics;

import etomica.api.IAction;
import etomica.api.IBox;
import etomica.api.IController;
import etomica.api.ISpecies;

import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.action.activity.Controller;
import etomica.modifier.ModifierNMolecule;
import etomica.simulation.prototypes.HSMD2D;

/**
 * Slider that selects the number of atoms of a given species in a box.
 *
 * @author David Kofke
 */
public class DeviceNSelector extends DeviceSlider {
    
    public DeviceNSelector() {
        this(null);
    }
    
    public DeviceNSelector(IController controller) {
        super(controller);
    }

    public void setResetAction(IAction newResetAction) {
        resetAction = newResetAction;
        if (modifyAction != null) {
            targetAction = new ActionGroupSeries(new IAction[]{modifyAction,resetAction});
        }
    }

    /**
     * Returns the action used to "reset" the simulation after changing the 
     * number of molecules, SimulationRestart by default.
     */
    public IAction getResetAction() {
        return resetAction;
    }

    public void setBox(IBox newBox) {
        box = newBox;
        if (species != null) {
            init();
        }
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
        if (box != null) {
            init();
        }
    }
    
    public IBox getBox() {
        return box;
    }
    
    public ISpecies getSpecies() {
        return species;
    }
    
    protected void init() {
        setMinimum(0);
        int max = 60;
        int nMolecules = box.getNMolecules(species);
        if (nMolecules > max) max = nMolecules;
        setMaximum(max);
        slider.setSnapToTicks(false);
        slider.setMajorTickSpacing(10);
        graphic(null).setSize(new java.awt.Dimension(40,30));
        setModifier(new ModifierNMolecule(box, species));
        if (resetAction != null) {
            targetAction = new ActionGroupSeries(new IAction[]{modifyAction,resetAction});
        }

        setLabel("Number of molecules");
    }
    
    protected IAction resetAction;
    protected ISpecies species;
    protected IBox box;
    
    //main method to demonstrate and test class
    public static void main(String[] args) {
        final String APP_NAME = "Devine n Selector";

        etomica.space.Space space = etomica.space2d.Space2D.getInstance();
        final HSMD2D sim = new HSMD2D();
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME, space);
        
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setBox(sim.box);
        nSelector.setSpecies(sim.species1);
        nSelector.setPostAction(graphic.getPaintAction(sim.box));
        graphic.add(nSelector);

        graphic.getController().getReinitButton().setPostAction(graphic.getPaintAction(sim.box));

        graphic.makeAndDisplayFrame(APP_NAME);

    }
    

} //end of DeviceNSelector
  