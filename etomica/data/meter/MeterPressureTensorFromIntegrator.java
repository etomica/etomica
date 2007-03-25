package etomica.data.meter;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.units.Pressure;

/**
 * Acts as a DataSource to retrieve the pressure from the integrator.
 * Currently only works with IntegratorVelocityVerlet (but could work just as
 * well with IntegratorVerlet;
 * see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=164 )
 */
public class MeterPressureTensorFromIntegrator implements DataSource, java.io.Serializable {

    public MeterPressureTensorFromIntegrator() {
        tag = new DataTag();
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setIntegrator(IntegratorVelocityVerlet newIntegrator) {
        integrator = newIntegrator;
        if (integrator != null) {
            data = new DataTensor(integrator.getPhase().getSpace());
            dataInfo = new DataTensor.DataInfoTensor("Pressure", Pressure.DIMENSION, integrator.getPhase().getSpace());
        }
    }
    
    public IntegratorPhase getIntegrator() {
        return integrator;
    }
    
    public Data getData() {
        data.x.E(integrator.getPressureTensor());
        return data;
    }
    
    private static final long serialVersionUID = 1L;
    protected IntegratorVelocityVerlet integrator;
    protected DataTensor data;
    protected DataInfo dataInfo;
    protected final DataTag tag;
}
