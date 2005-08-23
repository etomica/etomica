//This class includes a main method to demonstrate its use
package etomica.action;

import etomica.Action;
import etomica.EtomicaInfo;
import etomica.modifier.ModifierBoolean;

/**
 * Action that toggles a boolean value. This action can connect to any object
 * capable of switching between two states. The action operates through a
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object.
 * 
 * @author David Kofke
 */
public class ActionToggle implements Action, java.io.Serializable {

    private ModifierBoolean modifier;
    private String trueLabel = "True";
    private String falseLabel = "False";

    public ActionToggle(ModifierBoolean modifier) {
        this(modifier, "True", "False");
    }

    public ActionToggle(ModifierBoolean modifier, String trueText,
            String falseText) {
        setModifier(modifier);
        setTrueLabel(trueText);
        setFalseLabel(falseText);
    }

    /**
     * Toggles the state of the target object (true to false, false to true).
     */
    public void actionPerformed() {
        modifier.setBoolean(!modifier.getBoolean());
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Action that toggles a boolean value");
        return info;
    }

    /**
     * Returns the boolean modifier used by this action.
     */
    public ModifierBoolean getModifier() {
        return modifier;
    }

    /**
     * Specifies the boolean modifier that is set between true and false by the
     * action.
     */
    public void setModifier(ModifierBoolean newModifier) {
        modifier = newModifier;
    }

    /**
     * Specifies the button's label when the toggle is set to true.
     */
    public void setTrueLabel(String text) {
        trueLabel = text;
    }

    /**
     * @return the label displayed when the toggle is set to true.
     */
    public String getTrueLabel() {
        return trueLabel;
    }

    /**
     * Specifies the button's label when the toggle is set to false.
     */
    public void setFalseLabel(String text) {
        falseLabel = text;
    }

    /**
     * @return the label displayed when the toggle is set to false.
     */
    public String getFalseLabel() {
        return falseLabel;
    }

    /**
     * Returns the true or false label, depending on the current state of the
     * toggle.
     */
    public String getLabel() {
        return modifier.getBoolean() ? trueLabel : falseLabel;
    }
}
