package etomica.modules.statistics;

import etomica.data.AccumulatorAverage;
import etomica.data.IData;

public class DataAccSamplesPusher extends DataAccPusher {

    protected final AccumulatorAverage acc;

    public DataAccSamplesPusher(int idx, DataCollector collector, AccumulatorAverage acc) {
        super(idx, collector);
        this.acc = acc;
    }

    public void putData(IData data) {
        double stdev = data.getValue(0);
        double err = data.getValue(1);
        double nSamples = acc.getSampleCount() / ((stdev / err) * (stdev / err));
        collector.setData(idx, nSamples);
    }
}
