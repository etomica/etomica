package etomica.cavity;

import etomica.data.DataProcessorForked;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Null;

public class DataProcessorWidomContact extends DataProcessorForked {

    protected final DataGroup data;
    protected final DataDouble davg, derr;
    protected DataGroup contactData;

    public DataProcessorWidomContact() {
        super();
        DataDouble[] d = new DataDouble[3];
        DataDouble.DataInfoDouble[] di = new DataDouble.DataInfoDouble[3];
        for (int i = 0; i < 3; i++) {
            d[i] = new DataDouble();
            di[i] = new DataDouble.DataInfoDouble("mu bits", Null.DIMENSION);
        }
        davg = d[1];
        derr = d[2];
        d[0].x = Double.NaN;
        data = new DataGroup(d);
        dataInfo = new DataGroup.DataInfoGroup("mu", Null.DIMENSION, di);
    }

    @Override
    protected IData processData(IData inputData) {
        if (contactData == null) return null;
        DataGroup dg = (DataGroup) inputData;
        IData avg = dg.getData(0);
        IData err = dg.getData(1);
        double avgRatio = avg.getValue(avg.getLength() - 1);
        double errRatio = err.getValue(err.getLength() - 1);
        // contactData 0 is most recent
        double gcv = contactData.getData(1).getValue(0);
        double gce = contactData.getData(2).getValue(0);
        davg.x = Math.log(gcv / avgRatio);
        derr.x = Math.sqrt(gce * gce / (gcv * gcv) + errRatio * errRatio / (avgRatio * avgRatio));
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return dataInfo;
    }

    public GContactSink makeGContactSink() {
        return new GContactSink();
    }

    public class GContactSink implements IDataSink {

        @Override
        public void putData(IData data) {
            contactData = (DataGroup) data;
        }

        @Override
        public void putDataInfo(IDataInfo dataInfo) {
            if (!(dataInfo instanceof DataGroup.DataInfoGroup)) {
                throw new RuntimeException("oops");
            }
        }
    }
}
