package etomica;
import etomica.data.DataSourceAdapter;
import etomica.units.Dimension;

/**
 * A Modifier object permits changes to be made to a property of another object.
 * The property or object may be unknown to the object that is using the modifier.
 * For example, a modifier may be constructed to adjust a certain property, and then
 * the modifier can be passed to a generic Device, thereby connecting the device to the property.
 * @author David Kofke
 */

public abstract class Modifier extends DataSourceAdapter {
    
    public Modifier(Dimension dimension) {
        super(dimension);
    }
    
    public abstract void setValue(double d);
    public abstract double getValue();
    
    //to implement DataSource
    public double[] getData() {return new double[]{getValue()};}
    
}