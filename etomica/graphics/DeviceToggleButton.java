package etomica.graphics;

import java.awt.Color;

import etomica.Action;
import etomica.Controller;
import etomica.EtomicaInfo;
import etomica.action.ActionGroup;
import etomica.action.ActionToggle;
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

    public DeviceToggleButton(Controller controller, ModifierBoolean modifier) {
        this(controller, modifier, "True", "False");
    }

    public DeviceToggleButton(Controller controller, ModifierBoolean modifier,
            String trueText, String falseText) {
        super(controller);
        toggleAction = new ActionToggle(modifier, trueText, falseText);
        setAction(new ActionGroup(new Action[] { toggleAction, relabelButton }));
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Button that toggles a boolean value");
        return info;
    }

    /**
     * Toggles the button to the given state. No action is performed if the
     * button is already in the given state.
     */
    public void setState(boolean b) {
        if (b != getState())
            press();
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
            getButton().setText(getTrueLabel());
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
            getButton().setText(getFalseLabel());
    }

    /**
     * @return the label displayed when the toggle is set to false.
     */
    public String getFalseLabel() {
        return toggleAction.getFalseLabel();
    }

    private ActionToggle toggleAction;
    private final Action relabelButton = new Action() {

        public void actionPerformed() {
            DeviceToggleButton.this.getButton().setLabel(
                    toggleAction.getLabel());
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

        final etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);

        //here's the part unique to this class
        //sets up button to toggle atoms between red and blue
        ((DisplayPhase) graphic.displayList().getFirst())
                .setColorScheme(new ColorSchemeByType());
        ModifierBoolean modifier = new ModifierBoolean() {

            public void setBoolean(boolean b) {
                if (b)
                    ColorSchemeByType.setColor(sim.species, Color.RED);
                else
                    ColorSchemeByType.setColor(sim.species, Color.BLUE);
                //                sim.panel().repaint();
            }

            public boolean getBoolean() {
                return ColorSchemeByType.getColor(sim.species) == Color.RED;
            }
        };
        ActionToggle action = new ActionToggle(modifier);
        DeviceToggleButton button = new DeviceToggleButton(sim.getController(),
                modifier);
        button.setTrueLabel("Red");
        button.setFalseLabel("Blue");
        //end of unique part
        graphic.add(button);
        graphic.makeAndDisplayFrame();
    }

}