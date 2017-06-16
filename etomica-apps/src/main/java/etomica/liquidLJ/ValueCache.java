package etomica.liquidLJ;

import etomica.integrator.Integrator;
import etomica.data.DataSourceScalar;

public class ValueCache {
    protected long lastStep = -1;
    protected double lastValue;
    protected final DataSourceScalar dss;
    protected final Integrator integrator;
    public ValueCache(DataSourceScalar dss, Integrator integrator) {
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
