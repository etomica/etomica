package etomica;
import etomica.data.DataSourceAdapter;
import etomica.units.Dimension;

/**
 * A Modulator object permits changes to be made to a property of another object.
 * The property or object may be unknown to the object that is using the modulator.
 * For example, a modulator may be constructed to adjust a certain property, and then
 * the modulator can be passed to a generic Device, thereby connecting the device to the property.
 * @author David Kofke
 */

public abstract class ModulatorAbstract extends DataSourceAdapter {
    
    public ModulatorAbstract(Dimension dimension) {
        super(dimension);
    }
    
    public abstract void setValue(double d);
    public abstract double getValue();
    
    //to implement DataSource
    public double[] getData() {return new double[]{getValue()};}
    
}