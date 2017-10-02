/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;

import java.util.Arrays;


/**
 * Extracts one or more of the Data instances wrapped in a DataGroup.  If a single instance,
 * output is the instance itself; if multiple instances, a new DataGroup with the subset is
 * output.
 *
 * @author David Kofke
 *
 */

public class DataGroupFilter extends DataProcessor {

    private final boolean singleInstance;
    private final int[] indexes;
    private IData outputData;
    
    
    /**
     * Creates a filter that will take the DataGroup element indicated
     * by the index, counting from 0.
     */
    public DataGroupFilter(int index) {
        this(new int[] {index});
    }

    /**
     * Creates a filter that will take the set of elements from the DataGroup
     * corresponding to the given index values (with indexes counting from 0).
     */
    public DataGroupFilter(int[] indexes) {
        this.indexes = indexes.clone();
        singleInstance = (indexes.length == 1);
    }

    /**
     * Returns the output data that was established with a previous call to
     * processDataInfo.
     */
    public IData processData(IData data) {
        if (outputData == null) {
            if(singleInstance) {
                if(indexes[0] < ((DataGroup)data).getNData()) {
                    outputData = ((DataGroup)data).getData(indexes[0]);
                } else {
                    throw new ArrayIndexOutOfBoundsException("DataFilter was constructed to extract a Data element with an index that is larger than the number of Data elements wrapped in the DataGroup. Number of elements: "+((DataGroup)data).getNData()+"; index array: "+Arrays.toString(indexes));
                }
            } else {
                IData[] pushedData = new IData[indexes.length];
                try {
                    for (int i=0; i<indexes.length; i++) {
                        pushedData[i] = ((DataGroup)data).getData(indexes[i]);
                    }
                } catch(ArrayIndexOutOfBoundsException ex) {
                    throw new ArrayIndexOutOfBoundsException("DataFilter was constructed to extract a Data element with an index that is larger than the number of Data elements wrapped in the DataGroup. Number of elements: "+((DataGroup)data).getNData()+"; index array: "+Arrays.toString(indexes));
                }
                outputData = new DataGroup(pushedData);
            }
        }
        return outputData;
    }

    /**
     * Determines data that will be output, using the given DataInfo and
     * the index specification given at construction.  Returns the DataInfo
     * of the data that will be output.
     */
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        if (!(inputDataInfo instanceof DataInfoGroup)) {
            throw new IllegalArgumentException("DataGroupFilter must operate on a DataGroup");
        }
        outputData = null;
        int nData = ((DataInfoGroup)inputDataInfo).getNDataInfo();
        if(singleInstance) {
            if(indexes[0] < nData) {
                return ((DataInfoGroup)inputDataInfo).getSubDataInfo(indexes[0]);
            }
            throw new ArrayIndexOutOfBoundsException("DataFilter was constructed to extract a Data element with an index that is larger than the number of Data elements wrapped in the DataGroup. Number of elements: "+nData+"; index array: "+Arrays.toString(indexes));
        }
        IDataInfo[] pushedDataInfo = new IDataInfo[indexes.length];
        try {
            for (int i=0; i<indexes.length; i++) {
                pushedDataInfo[i] = ((DataInfoGroup)inputDataInfo).getSubDataInfo(indexes[i]);
            }
        } catch(ArrayIndexOutOfBoundsException ex) {
            throw new ArrayIndexOutOfBoundsException("DataFilter was constructed to extract a Data element with an index that is larger than the number of Data elements wrapped in the DataGroup. Number of elements: "+nData+"; index array: "+Arrays.toString(indexes));
        }
        DataInfoGroup myDataInfo = new DataInfoGroup(inputDataInfo.getLabel(), inputDataInfo.getDimension(), pushedDataInfo);
        myDataInfo.addTags(inputDataInfo.getTags());
        return myDataInfo;
    }
}
