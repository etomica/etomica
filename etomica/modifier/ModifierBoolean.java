package etomica.modifier;

/**
 * A ModifierBoolean object permits changes to be made to a boolean property of another object.
 */

public interface ModifierBoolean {
    
    public void setBoolean(boolean b);
    public boolean getBoolean();

}