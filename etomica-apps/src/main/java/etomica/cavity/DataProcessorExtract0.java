package etomica.cavity;

import etomica.data.DataProcessorForked;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Null;

class DataProcessorExtract0 extends DataProcessorForked {
    protected String label;
    protected final DataGroup data;
    protected DataProcessorFit errSource;

    public DataProcessorExtract0(String label, boolean hasErr) {
        this.label = label;
        DataDouble[] myData = new DataDouble[hasErr ? 3 : 1];
        myData[0] = new DataDouble();
        DataDouble.DataInfoDouble[] myDataInfo = new DataDouble.DataInfoDouble[myData.length];
        myDataInfo[0] = new DataDouble.DataInfoDouble("avg", Null.DIMENSION);
        if (myData.length == 3) {
            myData[1] = new DataDouble();
            myDataInfo[1] = new DataDouble.DataInfoDouble("the real avg", Null.DIMENSION);
            myData[2] = new DataDouble();
            myDataInfo[2] = new DataDouble.DataInfoDouble("err", Null.DIMENSION);
        }
        data = new DataGroup(myData);
        dataInfo = new DataGroup.DataInfoGroup(label, Null.DIMENSION, myDataInfo);
        dataInfo.addTag(tag);
    }

    public void setErrorSource(DataProcessorFit errSource) {
        this.errSource = errSource;
    }

    @Override
    protected IData processData(IData inputData) {
        if (errSource != null) {
            ((DataDouble) data.getData(0)).x = Double.NaN;
            ((DataDouble) data.getData(1)).x = inputData.getValue(0);
            ((DataDouble) data.getData(2)).x = errSource.getLastErr()[0];
        } else if (inputData instanceof DataGroup) {
            DataGroup dg = (DataGroup) inputData;
            ((DataDouble) data.getData(0)).x = dg.getData(0).getValue(0);
            ((DataDouble) data.getData(1)).x = dg.getData(1).getValue(0);
            ((DataDouble) data.getData(2)).x = dg.getData(2).getValue(0);
        } else {
            ((DataDouble) data.getData(0)).x = inputData.getValue(0);
        }
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return dataInfo;
    }
}
