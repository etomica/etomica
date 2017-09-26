package etomica.modules.statistics;

import etomica.data.AccumulatorAverage;
import etomica.data.IData;

public class DataAccDifficultyPusher extends DataAccPusher {

    protected final AccumulatorAverage acc;

    public DataAccDifficultyPusher(int idx, DataCollector collector, AccumulatorAverage acc) {
        super(idx, collector);
        this.acc = acc;
    }

    public void putData(IData data) {
        collector.setData(idx, data.getValue(0) * Math.sqrt(acc.getSampleCount()));
    }
}
