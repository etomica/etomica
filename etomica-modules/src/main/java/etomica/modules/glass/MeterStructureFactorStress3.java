package etomica.modules.glass;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.meter.MeterStructureFactor;
import etomica.data.types.DataFunction;

public class MeterStructureFactorStress3 implements IDataSource {
    private final MeterStructureFactor[] meters;
    private final DataTag tag;
    private final IDataInfo dataInfo;
    private final DataFunction data;

    public MeterStructureFactorStress3(MeterStructureFactor[] meters) {
        this.meters = meters;
        tag = new DataTag();
        dataInfo = meters[0].getDataInfo().getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = (DataFunction) dataInfo.makeData();
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        double[] y = data.getData();
        IData[] meterData = new IData[meters.length];
        for (int j = 0; j < meters.length; j++) {
            meterData[j] = meters[j].getData();
        }
        for (int i = 0; i < y.length; i++) {
            y[i] = 0;
            for (int j = 0; j < meters.length; j++) {
                y[i] += meterData[j].getValue(i);
            }
        }
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
}
