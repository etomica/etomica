package etomica;

/**
 * A ModulatorBoolean object permits changes to be made to a boolean property of another object.
 */

public abstract class ModulatorBoolean implements java.io.Serializable {
    
    public abstract void setBoolean(boolean b);
    public abstract boolean getBoolean();

}