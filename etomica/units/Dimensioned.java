package etomica.units;

/**
 * Interface for an object with an associated physical dimension.
 * Indicates that the object has methods for setting and getting its units.
 */

public interface Dimensioned {
    
    public void setUnit(Unit u);
    public Unit getUnit();
    public Dimension dimension();
    
}