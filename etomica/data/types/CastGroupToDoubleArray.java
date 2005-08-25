package etomica.data.types;


import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataGroupExtractor;
import etomica.data.DataInfo;
import etomica.data.DataJudge;
import etomica.data.DataProcessor;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * A DataProcessor that converts a homogeneous DataGroup into a multidimensional DataDoubleArray. 
 * All elements in group must be of the same data type.  If the group contains only one element,
 * the effect is the same as applying CastToDoubleArray to the element.  Otherwise, if <tt>d</tt>
 * is the dimension of each grouped element, the resulting DataDoubleArray will be of 
 * dimension <tt>d+1</tt>.  The added dimension corresponds to the different elements of the group.
 * <p>
 * For example:
 * <ul>
 * <li>if the DataGroup contains 8 DataDouble instances (thus d = 0), the cast gives a
 * one-dimensional DataDoubleArray of length 8.
 * <li>if the DataGroup contains 4 Vector3D instances (thus d = 1), the cast gives a two-dimensional DataDoubleArray
 * of dimensions (4, 3).
 * <li>if the DataGroup contains 5 two-dimensional (d = 2) DataDoubleArray of shape (7,80), the cast gives a
 * three-dimensional DataDoubleArray of shape (5, 7, 80).
 * </ul>
 * etc. 
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastGroupToDoubleArray extends DataProcessor {

    /**
     * Sole constructor.
     */
    public CastGroupToDoubleArray() {
    }

    /**
     * Prepares processor to handle Data. Uses given DataInfo to determine the
     * type of Data to expect in subsequent calls to processData.
     * 
     * @throws IllegalArgumentException
     *             if DataInfo does not indicate a DataGroup, or if DataInfo
     *             indicates that expected DataGroup will not be homogeneous
     */
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
                outputData = new DataDoubleArray(label, dimension, size);
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
     * Converts data in given group to a DataDoubleArray as described in the
     * general comments for this class.
     * 
     * @throws ClassCastException
     *             if the given Data is not a DataGroup with Data elements of
     *             the type indicated by the most recent call to
     *             processDataInfo.
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
    
    /**
     * Returns null.
     */
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataDoubleArray outputData;
    private int inputType;
}
