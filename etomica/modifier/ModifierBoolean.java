package etomica.modifier;

/**
 * A ModifierBoolean object permits changes to be made to a boolean property of another object.
 */

public abstract class ModifierBoolean implements java.io.Serializable {
    
    public abstract void setBoolean(boolean b);
    public abstract boolean getBoolean();

}