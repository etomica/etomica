package etomica.data.types;


import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.lattice.IndexIteratorSequential;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * A DataProcessor that converts a DataArray instance into a DataDoubleArray. 
 * Copies an element of the input data to the DataDouble's encapsulated value 
 * and returns the DataDoubleArray instance.  No attempt is made to reduce the 
 * number of dimenions in the resulting DataDoubleArray.  If the incoming 
 * DataArray is 1x1x1 and wraps a DataDouble, the resulting DataDoubleArray 
 * also be 3 dimensional.
 * <p> 
 * Can cast from DataDouble, DataDoubleArray, DataInteger, DataVector, DataTensor
 * and DataFunction classes wrapped in a DataArray.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastArrayToDoubleArray extends DataProcessor {

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if (inputDataInfo.getDataClass() != DataArray.class) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        D = 0;
        String label = inputDataInfo.getLabel();
        Dimension dimension = inputDataInfo.getDimension();
        DataArray.Factory factory = (DataArray.Factory)inputDataInfo.getDataFactory();
        DataFactory subFactory = factory.getArrayElementFactory();
//        Data[] data = ((DataArray.Factory)inputDataInfo.getDataFactory()).data;
        if (factory.getArrayLength() == 0) {
            inputType = 0;
            outputData = new DataDoubleArray(label,dimension,0);
            return outputData.getDataInfo();
        }
        Class innerDataClass = subFactory.getDataClass();
        int[] outerArraySize = factory.getArraySize();
        indexIterator = new IndexIteratorSequential(outerArraySize.length);
        indexIterator.setSize(outerArraySize);
        
        if (innerDataClass == DataDoubleArray.class) {
            inputType = 1;
            int[] size = ((DataDoubleArray.Factory)subFactory).getArrayShape();
            int[] totalSize = new int[size.length+outerArraySize.length];
            System.arraycopy(outerArraySize,0,totalSize,0,outerArraySize.length);
            System.arraycopy(size,0,totalSize,outerArraySize.length,size.length);
            outputData = new DataDoubleArray(label,dimension,totalSize);
        } else if (innerDataClass == DataDouble.class) {
            inputType = 2;
            outputData = new DataDoubleArray(label, dimension, outerArraySize);
        } else if (innerDataClass == DataInteger.class) {
            inputType = 3;
            outputData = new DataDoubleArray(label, dimension, outerArraySize);
        } else if (innerDataClass == DataVector.class) {
            inputType = 4;
            int[] totalSize = new int[outerArraySize.length+1];
            System.arraycopy(outerArraySize,0,totalSize,0,outerArraySize.length);
            D = ((DataVector.Factory)subFactory).getSpace().D();
            totalSize[totalSize.length-1] = D;
            outputData = new DataDoubleArray(label,dimension,totalSize);
        } else if (innerDataClass == DataTensor.class) {
            inputType = 5;
            int[] totalSize = new int[outerArraySize.length+2];
            System.arraycopy(outerArraySize,0,totalSize,0,outerArraySize.length);
            D = ((DataTensor.Factory)subFactory).getSpace().D();
            totalSize[totalSize.length-2] = D;
            totalSize[totalSize.length-1] = D;
            outputData = new DataDoubleArray(label,dimension,totalSize);
        } else if(innerDataClass == DataFunction.class) {
            inputType = 6;
            int[] sizes =  ((DataFunction.Factory)subFactory).getIndependentDataSizes();
            
            int[] totalSize = new int[sizes.length+outerArraySize.length];
            System.arraycopy(outerArraySize,0,totalSize,0,outerArraySize.length);
            System.arraycopy(sizes,0,totalSize,outerArraySize.length,sizes.length);
            outputData = new DataDoubleArray(label,dimension,totalSize);
        } else if(innerDataClass == DataGroup.class) {
            throw new IllegalArgumentException("Ok, now you're just being evil; I don't know how to cast a DataArray-wrapped DataGroup into a DataDoubleArray.  Unwrap your groups before sticking them in a DataArray perhaps");
        } else {
            throw new IllegalArgumentException("Cannot cast to DataDoubleArray from "
                    + innerDataClass + " in DataArray");
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
        DataArray dataArray= (DataArray)data;
        indexIterator.reset();
        double[] x = outputData.getData();
        int i = 0;
        switch (inputType) {
        case 0:  // empty group 
            return outputData;
        case 1:  // DataDoubleArray
            while (indexIterator.hasNext()) {
                int[] idx = indexIterator.next();
                outputData.assignSubsectionFrom(idx,((DataDoubleArray)dataArray.getData(idx)).getData());
            }
            return outputData;
        case 2:  // DataDouble(s)
            while (indexIterator.hasNext()) {
                int[] idx = indexIterator.next();
                x[i++] = ((DataDouble)dataArray.getData(idx)).x;
            }
            return outputData;
        case 3:  // DataInteger(s)
            while (indexIterator.hasNext()) {
                int[] idx = indexIterator.next();
                x[i++] = ((DataDouble)dataArray.getData(idx)).x;
            }
            return outputData;
        case 4:  // a single DataVector
            while (indexIterator.hasNext()) {
                int[] idx = indexIterator.next();
                Vector v = ((DataVector)dataArray.getData(idx)).x;
                for (int j=0; j<D; j++) {
                    x[i++] = v.x(j);
                }
            }
            return outputData;
        case 5:  // a single DataTensor
            while (indexIterator.hasNext()) {
                int[] idx = indexIterator.next();
                Tensor t = ((DataTensor)dataArray.getData(idx)).x;
                for (int j=0; j<D; j++) {
                    for (int l=0; l<D; l++) {
                        x[i++] = t.component(j,l);
                    }
                }
            }
            return outputData;
        case 6:
            while (indexIterator.hasNext()) {
                int[] idx = indexIterator.next();
                outputData.assignSubsectionFrom(idx,((DataFunction)dataArray.getData(idx)).getYData().getData());
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

    private IndexIteratorSequential indexIterator;
    private DataDoubleArray outputData;
    private int inputType;
    private int D;
}
