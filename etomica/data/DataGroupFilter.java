package etomica.data;

import etomica.Data;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

public class DataGroupFilter extends DataPipe {

    public DataGroupFilter(int[] indexes) {
        this.indexes = (int[])indexes.clone();
    }

    public void putData(Data data) {
        DataGroup dataGroup = (DataGroup)data;
        if (pushedDataGroup == null) {
            initialize(dataGroup, indexes);
        }
        for(int i=0; i<indexes.length; i++) {
            pushedDataGroup.getData(i).E(dataGroup.getData(indexes[i]));
        }
        pushData(pushedDataGroup);
    }
    
    private void initialize(DataGroup dataGroup, int[] indexes) {
        Data[] pushedData = new Data[indexes.length];
        for (int i=0; i<indexes.length; i++) {
            pushedData[i] = dataGroup.getData(indexes[i]).makeCopy();
        }
        pushedDataGroup = new DataGroup(dataGroup.getDataInfo().getLabel(),
                                        dataGroup.getDataInfo().getDimension(),
                                        pushedData);
    }
    
    private DataGroup pushedDataGroup;
    private final int[] indexes;
}
