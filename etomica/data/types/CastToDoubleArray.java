package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataFactory;
import etomica.data.DataGroupExtractor;
import etomica.data.DataJudge;
import etomica.data.DataProcessor;
import etomica.units.Dimension;

/**
 * A DataProcessor that converts a Data instance into a DataDoubleArray. Copies an
 * element of the input data to the DataDouble's encapsulated value and returns
 * the DataDouble instance.
 * <p> 
 * Can cast from DataDouble, DataDoubleArray, DataInteger, DataVector, DataTensor.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastToDoubleArray extends DataProcessor {

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        Class inputClass = inputDataInfo.getClass();
        String label = inputDataInfo.getLabel();
        Dimension dimension = inputDataInfo.getDimension();
        DataFactory factory = inputDataInfo.getDataFactory();
        if (inputClass == DataDoubleArray.class) {
            inputType = 0;
            return inputDataInfo;
        } else if (inputClass == DataDouble.class) {
            inputType = 1;
            outputData = new DataDoubleArray(label, dimension, 1);
        } else if (inputClass == DataInteger.class) {
            inputType = 2;
            outputData = new DataDoubleArray(label, dimension, 1);
        } else if (inputClass == DataVector.class) {
            inputType = 3;
            outputData = new DataDoubleArray(label, dimension, ((DataVector.Factory)factory).getSpace().D());
        } else if (inputClass == DataTensor.class) {
            inputType = 4;
            int D = ((DataTensor.Factory)factory).getSpace().D();
            outputData = new DataDoubleArray(label, dimension, new int[] {D, D}); 
        } else if(inputClass == DataFunction.class) {
            inputType = 5;
            outputData = new DataDoubleArray(label, dimension, ((DataFunction.Factory)factory).getIndependentDataSizes());
        } else if(inputClass == DataGroup.class) {
            inputType = 6;
            DataGroupExtractor extractor = new DataGroupExtractor(new DataJudge.ByClass(DataDoubleArray.class, true));
            extractor.putDataInfo(inputDataInfo);
            outputData = (DataDoubleArray)extractor.processData(null);
        } else {
            throw new IllegalArgumentException("Cannot cast to DataDoubleArray from "
                    + inputClass);
        }
        return outputData.getDataInfo();
    }
    
    /**
     * Extracts a double from the input data and returns it encapsulated in a
     * DataDouble.
     * 
     * @param data
     *            a Data instance of the type indicated by the DataInfo at
     *            construction
     * @return a DataDouble holding the value cast from the given Data; the same
     *         instance is returned with every invocation.
     * 
     * @throws ClassCastException
     *             if the given Data is not of the same type as indicated by the
     *             DataInfo given at construction
     */
    protected Data processData(Data data) {
        switch (inputType) {
        case 0:
            return data;
        case 1:
            outputData.E(((DataDouble) data).x);
            return outputData;
        case 2:
            outputData.E(((DataInteger) data).x);
            return outputData;
        case 3:
            ((DataVector)data).x.assignTo(outputData.getData());
            return outputData;
        case 4:
            ((DataTensor)data).x.assignTo(outputData.getData());//both Tensor and DataDoubleArray sequence data by rows
            return outputData;
        case 5:
            outputData.E(((DataFunction)data).getYData());
            return outputData;
        case 6:
            return outputData;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataDoubleArray outputData;
    private int inputType;
}
