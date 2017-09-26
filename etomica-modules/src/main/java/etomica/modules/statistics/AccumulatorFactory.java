package etomica.modules.statistics;

import etomica.data.AccumulatorAverageFixed;

public interface AccumulatorFactory {
    public AccumulatorAverageFixed makeAccumulator();
}
