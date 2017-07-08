/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;


import etomica.space.Vector;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.data.types.DataVector.DataInfoVector;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;

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
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        if (!(inputDataInfo instanceof DataInfoGroup)) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        String label = inputDataInfo.getLabel();
        Dimension dimension = inputDataInfo.getDimension();
        int numSubData = ((DataInfoGroup)inputDataInfo).getNDataInfo();
        if (numSubData == 0) {
            inputType = 0;
            outputData = new DataDoubleArray(0);
            IEtomicaDataInfo outputDataInfo = new DataInfoDoubleArray(label, dimension, new int[]{0});
            outputDataInfo.addTags(inputDataInfo.getTags());
            return outputDataInfo;
        }
        IEtomicaDataInfo subDataInfo = ((DataInfoGroup)inputDataInfo).getSubDataInfo(0);

        Class subDataInfoClass = subDataInfo.getClass();
        for (int i = 1; i<numSubData; i++) {
            if (((DataInfoGroup)inputDataInfo).getSubDataInfo(i).getClass() != subDataInfoClass) {
                throw new IllegalArgumentException("CastGroupToDoubleArray can only handle homogeneous groups");
            }
        }
        
        int[] outputArrayShape;
        if (subDataInfo instanceof DataInfoDoubleArray) {
            int[] arrayShape = ((DataInfoDoubleArray)subDataInfo).getArrayShape();
            int D = arrayShape.length;
            if (numSubData == 1) {
                outputArrayShape = arrayShape;
                inputType = 1;
                // we don't need an outputData instance
                return subDataInfo;
            }
            outputArrayShape = new int[++D];
            outputArrayShape[0] = numSubData;
            for (int i=1; i<D; i++) {
                outputArrayShape[i] = arrayShape[i-1];
            }
            inputType = 2;
        } else if (subDataInfo instanceof DataInfoDouble) {
            inputType = 3;
            outputArrayShape = new int[]{numSubData};
        } else if (subDataInfo instanceof DataInfoVector) {
            int D = ((DataInfoVector)subDataInfo).getSpace().D();
            if (numSubData == 1) {
                inputType = 4;
                outputArrayShape = new int[]{D};
            }
            else {
                inputType = 5;
                outputArrayShape = new int[]{numSubData, D};
            }
        } else if (subDataInfo instanceof DataInfoTensor) {
            int D = ((DataInfoTensor)subDataInfo).getSpace().D();
            if (numSubData == 1) {
                inputType = 6;
                outputArrayShape = new int[]{D,D};
            }
            else {
                inputType = 7;
                outputArrayShape = new int[]{numSubData,D,D};
            }
        } else {
            throw new IllegalArgumentException("Cannot cast to DataDoubleArray from "
                    + subDataInfo.getClass() + " in DataGroup");
        }
        outputData = new DataDoubleArray(outputArrayShape);
        DataInfoDoubleArray outputDataInfo = new DataInfoDoubleArray(label, dimension, outputArrayShape);
        outputDataInfo.addTags(inputDataInfo.getTags());
        return outputDataInfo;
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
    protected IData processData(IData data) {
        DataGroup group = (DataGroup)data;
        //we don't add ourselves
        switch (inputType) {
        case 0:  // empty group 
            return outputData;
        case 1:  // a single DataDoubleArray
            return ((DataGroup)data).getData(0);
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
        case 4:  // a single DataVector
            group.getData(0).assignTo(outputData.getData());
            return outputData;
        case 5:  // multiple DataVectors
            double[] x = outputData.getData();
            int k=0;
            for (int i=0; i<group.getNData(); i++) {
                Vector v = ((DataVector)group.getData(i)).x;
                for (int j=0; j<v.getD(); j++) {
                    x[k++] = v.getX(j);
                }
            }
            return outputData;
        case 6:  // a single DataTensor
            group.getData(0).assignTo(outputData.getData());
            return outputData;
        case 7:  // multiple DataTensors
            x = outputData.getData();
            k=0;
            for (int i=0; i<group.getNData(); i++) {
                Tensor t = ((DataTensor)group.getData(i)).x;
                for (int j=0; j<t.D(); j++) {
                    for (int l=0; l<t.D(); l++) {
                        x[k++] = t.component(j,l);
                    }
                }
            }
            return outputData;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    /**
     * Returns null.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo info) {
        return null;
    }

    private static final long serialVersionUID = 1L;
    private DataDoubleArray outputData;
    private int inputType;
}
