package etomica.modules.adsorption;

import etomica.api.IBoundary;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.ISpace;

public class MeterProfileByVolumeAdsorption extends MeterProfileByVolume {

    public MeterProfileByVolumeAdsorption(ISpace space) {
        super(space);
        setProfileDim(1);
    }
    
    public void reset() {
        if (box == null || xMax-xMin==0) return;
        
        IBoundary boundary = box.getBoundary();
        xDataSource.setXMin(xMin);
        xDataSource.setXMax(xMax);
        xDataSource.setDoEnforceBounds(false);
        
        if (meter != null) {
            data = new DataFunction(new int[] {xDataSource.getNValues()});
            dataInfo = new DataInfoFunction(meter.getMoleculeDataInfo().getLabel()+" Profile", meter.getMoleculeDataInfo().getDimension(), this);
            dataInfo.addTag(meter.getTag());
            dataInfo.addTag(tag);

            dV = (xMax-xMin)/data.getLength();
            for (int i=0; i<boundary.getBoxSize().getD(); i++) {
                if (i==profileDim) continue;
                dV *= boundary.getBoxSize().getX(i);
            }
        }
    }

    public void setRange(double xMin, double xMax) {
        this.xMin = xMin;
        this.xMax = xMax;
        reset();
    }

    protected double xMin, xMax;
}
