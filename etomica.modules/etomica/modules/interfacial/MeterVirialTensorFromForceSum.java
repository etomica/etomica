package etomica.modules.interfacial;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.space.ISpace;
import etomica.units.Energy;

/**
 * Acts as a DataSource to retrieve the pressure from the integrator.
 * Currently only works with IntegratorVelocityVerlet (but could work just as
 * well with IntegratorVerlet;
 * see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=164 )
 */
public class MeterVirialTensorFromForceSum implements DataSource, java.io.Serializable {

    public MeterVirialTensorFromForceSum(ISpace space, PotentialCalculationForcePressureBinSum forceSum) {
        tag = new DataTag();
        this.forceSum = forceSum;
        data = new DataTensor(space);
        dataInfo = new DataTensor.DataInfoTensor("Virial", Energy.DIMENSION, space);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        data.x.E(forceSum.getPressureTensor());
        return data;
    }
    
    private static final long serialVersionUID = 1L;
    protected PotentialCalculationForcePressureBinSum forceSum;
    protected DataTensor data;
    protected DataInfo dataInfo;
    protected final DataTag tag;
}
