package etomica.modifier;

/**
 * Interface that permits changes to be made to a boolean property of another object.
 * This type of modifier might be plugged into a checkbox or toggle button.
 */

public interface ModifierBoolean {
    
    /**
     * Sets the modified value to the given boolean.
     */
    public void setBoolean(boolean b);

    /**
     * Returns the current state of the modified value.
     */
    public boolean getBoolean();

}