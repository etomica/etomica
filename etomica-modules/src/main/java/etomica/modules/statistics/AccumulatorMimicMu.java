package etomica.modules.statistics;

import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;

public class AccumulatorMimicMu extends AccumulatorAverageFixed {

    protected final IntegratorBox integrator;

    public AccumulatorMimicMu(IntegratorBox integrator) {
        super();
        this.integrator = integrator;
    }

    public void putDataInfo(IDataInfo dataInfo) {
        dataSink.putDataInfo(dataInfo);
    }

    public IData processData(IData data) {
        double temperature = integrator.getTemperature();
        DataGroup g = (DataGroup) data;
        DataDouble avg = (DataDouble) g.getData(AVERAGE.index);
        DataDouble err = (DataDouble) g.getData(ERROR.index);
        DataDouble stdev = (DataDouble) g.getData(STANDARD_DEVIATION.index);
        err.x = temperature * err.x / avg.x;
        stdev.x = temperature * stdev.x / avg.x;
        avg.x = -temperature * Math.log(avg.x);
        return data;
    }
}
