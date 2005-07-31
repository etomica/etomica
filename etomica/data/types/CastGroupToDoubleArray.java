package etomica.data.types;


import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataFactory;
import etomica.data.DataGroupExtractor;
import etomica.data.DataJudge;
import etomica.data.DataProcessor;
import etomica.space.Tensor;
import etomica.space.Vector;
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
public class CastGroupToDoubleArray extends DataProcessor {

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if (inputDataInfo.getDataClass() != DataGroup.class) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        String label = inputDataInfo.getLabel();
        Dimension dimension = inputDataInfo.getDimension();
        DataFactory factory = inputDataInfo.getDataFactory();
        Data[] data = ((DataGroup.Factory)inputDataInfo.getDataFactory()).data;
        if (data.length == 0) {
            inputType = 0;
            outputData = new DataDoubleArray(label,dimension,0);
            return outputData.getDataInfo();
        }
        Class innerDataClass = data[0].getClass();
        for (int i = 1; i<data.length; i++) {
            if (data[i].getClass() != innerDataClass) {
                throw new IllegalArgumentException("CastGroupToDoubleArray can only handle homogeneous groups");
            }
        }
        
        if (innerDataClass == DataDoubleArray.class) {
            DataDoubleArray dataArray = (DataDoubleArray)data[0];
            int D = dataArray.getArrayDimension();
            int[] size;
            if (data.length == 1) {
                inputType = 1;
                outputData = dataArray;
            }
            else {
                size = new int[++D];
                size[0] = data.length;
                for (int i=1; i<D; i++) {
                    size[i] = dataArray.getArrayShape(i-1);
                }
                inputType = 2;
            }
        } else if (innerDataClass == DataDouble.class) {
            inputType = 3;
            outputData = new DataDoubleArray(label, dimension, data.length);
        } else if (innerDataClass == DataInteger.class) {
            inputType = 4;
            outputData = new DataDoubleArray(label, dimension, data.length);
        } else if (innerDataClass == DataVector.class) {
            int D = ((DataVector.Factory)factory).getSpace().D();
            if (data.length == 1) {
                inputType = 5;
                outputData = new DataDoubleArray(label, dimension, D);
            }
            else {
                inputType = 6;
                outputData = new DataDoubleArray(label, dimension, new int[]{data.length,D});
            }
        } else if (innerDataClass == DataTensor.class) {
            int D = ((DataTensor.Factory)factory).getSpace().D();
            if (data.length == 1) {
                inputType = 7;
                outputData = new DataDoubleArray(label, dimension, new int[] {D, D});
            }
            else {
                inputType = 8;
                outputData = new DataDoubleArray(label, dimension, new int[] {data.length, D, D});
            }
        } else if(innerDataClass == DataFunction.class) {
            int[] sizes =  ((DataFunction.Factory)factory).getIndependentDataSizes();
            inputType = 9;
            if (data.length != 1) {
                inputType = 10;
                int[] newSizes = new int[sizes.length+1];
                newSizes[0] = data.length;
                System.arraycopy(sizes,0,newSizes,1,sizes.length);
                sizes = newSizes;
            }
            outputData = new DataDoubleArray(label, dimension, sizes);
        } else if(innerDataClass == DataGroup.class) {
            inputType = 11;
            DataGroupExtractor extractor = new DataGroupExtractor(new DataJudge.ByClass(DataDoubleArray.class, true));
            extractor.putDataInfo(inputDataInfo);
            outputData = (DataDoubleArray)extractor.processData(null);
        } else {
            throw new IllegalArgumentException("Cannot cast to DataDoubleArray from "
                    + innerDataClass + "in DataGroup");
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
        DataGroup group = (DataGroup)data;
        switch (inputType) {
        case 0:  // empty group 
            return outputData;
        case 1:  // a single DataDoubleArray
            return outputData;
        case 2:  // multiple DataDoubleArrays
            for (int i=0; i<group.getNData(); i++) {
                outputData.assignColumnFrom(i,((DataDoubleArray)group.getData(i)).getData());
            }
            return outputData;
        case 3:  // DataDouble(s)
            for (int i=0; i<group.getNData(); i++) {
                outputData.getData()[i] = ((DataDouble)group.getData(i)).x;
            }
            return outputData;
        case 4:  // DataInteger(s)
            for (int i=0; i<group.getNData(); i++) {
                outputData.getData()[i] = ((DataInteger)group.getData(i)).x;
            }
            return outputData;
        case 5:  // a single DataVector
            ((DataVector)group.getData(0)).assignTo(outputData.getData());
            return outputData;
        case 6:  // multiple DataVectors
            double[] x = outputData.getData();
            int k=0;
            for (int i=0; i<group.getNData(); i++) {
                Vector v = ((DataVector)group.getData(i)).x;
                for (int j=0; j<v.D(); j++) {
                    x[k++] = v.x(j);
                }
            }
            return outputData;
        case 7:  // a single DataTensor
            ((DataTensor)group.getData(0)).assignTo(outputData.getData());
            return outputData;
        case 8:  // multiple DataTensors
            x = outputData.getData();
            k=0;
            for (int i=0; i<group.getNData(); i++) {
                Tensor t = ((DataTensor)group.getData(i)).x;
                for (int j=0; j<t.length(); j++) {
                    for (int l=0; l<t.length(); l++) {
                        x[k++] = t.component(j,l);
                    }
                }
            }
            return outputData;
        case 9:
            outputData.E(((DataFunction)group.getData(0)).getYData());
            return outputData;
        case 10:
            for (int i=0; i<group.getNData(); i++) {
                outputData.assignColumnFrom(i,((DataFunction)group.getData(0)).getYData().getData());
            }
            return outputData;
        case 11:
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
