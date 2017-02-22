package etomica.liquidLJ;

import etomica.api.IIntegrator;
import etomica.data.DataSourceScalar;

public class ValueCache {
    protected long lastStep = -1;
    protected double lastValue;
    protected final DataSourceScalar dss;
    protected final IIntegrator integrator;
    public ValueCache(DataSourceScalar dss, IIntegrator integrator) {
        this.dss = dss;
        this.integrator = integrator;
    }
    public double getValue() {
        if (integrator.getStepCount() != lastStep) {
            lastStep = integrator.getStepCount();
            lastValue = dss.getDataAsScalar();
        }
        return lastValue;
    }
}