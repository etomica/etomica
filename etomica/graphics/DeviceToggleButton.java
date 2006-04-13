package etomica.graphics;

import java.awt.Color;

import etomica.EtomicaInfo;
import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.ActionToggle;
import etomica.action.activity.Controller;
import etomica.modifier.ModifierBoolean;

/**
 * Button that toggles a boolean value. This device can connect to any object
 * capable of switching between two states. The device operates through a
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object. Action is performed through the controller if given at
 * construction or set afterwards.
 * 
 * @author David Kofke
 */
public class DeviceToggleButton extends DeviceButton {

    public DeviceToggleButton(Controller controller) {
        super(controller);
    }
    
    public DeviceToggleButton(Controller controller, ModifierBoolean modifier) {
        this(controller, modifier, "True", "False");
    }

    public DeviceToggleButton(Controller controller, ModifierBoolean modifier,
            String trueText, String falseText) {
        super(controller);
        setModifier(modifier, trueText, falseText);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Button that toggles a boolean value");
        return info;
    }
    
    public void setModifier(ModifierBoolean modifier, String trueText, String falseText) {
        toggleAction = new ActionToggle(modifier, trueText, falseText);
        setAction(new ActionGroupSeries(new Action[] { toggleAction, relabelButton }));
        button.setText(modifier.getBoolean() ? trueText : falseText);
    }

    /**
     * Toggles the button to the given state. No action is performed if the
     * button is already in the given state.
     */
    public void setState(boolean b) {
        if (b != getState())
            button.doClick();
    }

    /**
     * Returns the current true/false state of the button.
     */
    public boolean getState() {
        return toggleAction.getModifier().getBoolean();
    }

    /**
     * Specifies the button's label when the toggle is set to true.
     */
    public void setTrueLabel(String text) {
        toggleAction.setTrueLabel(text);
        if (getState())
            button.setText(getTrueLabel());
    }

    /**
     * @return the label displayed when the toggle is set to true.
     */
    public String getTrueLabel() {
        return toggleAction.getTrueLabel();
    }

    /**
     * Specifies the button's label when the toggle is set to false.
     */
    public void setFalseLabel(String text) {
        toggleAction.setFalseLabel(text);
        if (!getState())
            button.setText(getFalseLabel());
    }

    /**
     * @return the label displayed when the toggle is set to false.
     */
    public String getFalseLabel() {
        return toggleAction.getFalseLabel();
    }

    protected ActionToggle toggleAction;
    private final Action relabelButton = new Action() {

        public void actionPerformed() {
            DeviceToggleButton.this.setLabel(toggleAction.getLabel());
        }

        public String getLabel() {
            return "";
        }
    };

    /**
     * Method to demonstrate and test the use of this class. Slider is used to
     * control the temperature of a hard-sphere MD simulation
     */
    public static void main(String[] args) {

        final etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);

        //here's the part unique to this class
        //sets up button to toggle atoms between red and blue
        final ColorSchemeByType colorScheme = new ColorSchemeByType();
        ((DisplayPhase) graphic.displayList().getFirst())
                .setColorScheme(colorScheme);
        ModifierBoolean modifier = new ModifierBoolean() {

            public void setBoolean(boolean b) {
                if (b)
                    colorScheme.setColor(sim.species.getMoleculeType(), Color.RED);
                else
                    colorScheme.setColor(sim.species.getMoleculeType(), Color.BLUE);
                //                sim.panel().repaint();
            }

            public boolean getBoolean() {
                return colorScheme.getColor(sim.species.getMoleculeType()) == Color.RED;
            }
        };
        DeviceToggleButton button = new DeviceToggleButton(sim.getController(),
                modifier);
        button.setTrueLabel("Red");
        button.setFalseLabel("Blue");
        //end of unique part
        graphic.add(button);
        graphic.makeAndDisplayFrame();
    }

}