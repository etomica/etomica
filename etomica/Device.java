package etomica;
import etomica.units.Unit;
import etomica.units.Dimensioned;
import etomica.units.Dimension;
import java.awt.Component;

/**
 * Base class for all Devices.  These are objects that permit manipulation of the
 * fields of other objects.  They provide a (usually) graphical interface to permit
 * values of fields to be input and changed.
 */
public abstract class Device extends SimulationElement implements Simulation.GraphicalElement, Dimensioned, java.io.Serializable {
    
    public static final String VERSION = "Device:01.01.17";
    
    protected Unit unit;
    
    public Device(Simulation sim) {
        super(sim, Device.class);
    }
        
    public abstract Component graphic(Object obj);
    
    public void setUnit(Unit u) {unit = u;}
    public Unit getUnit() {return unit;}
    public Dimension dimension() {return unit.dimension();} //may want to override this in most cases
    
}