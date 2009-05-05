package etomica.data.meter;

import etomica.api.IBox;
import etomica.data.DataSourceScalar;
import etomica.units.Volume;

/**
 * Meter for measurement of the box volume.
 */
public class MeterVolume extends DataSourceScalar {
    
    public MeterVolume() {
        super("Volume",Volume.DIMENSION);
    }

    public double getDataAsScalar() {
        return box.getBoundary().volume();
    }
    
    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
}
