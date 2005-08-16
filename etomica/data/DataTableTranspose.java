package etomica.data;

import etomica.data.types.CastToTable;
import etomica.data.types.DataTable;


/**
 * Takes a DataTable and returns a different DataTable obtained by performing
 * a matrix transpose of the input table elements.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jul 28, 2005 by kofke
 */
public class DataTableTranspose extends DataProcessor {

    /* (non-Javadoc)
     * @see etomica.data.DataProcessor#processData(etomica.Data)
     */
    protected Data processData(Data inputData) {
        DataTable inputTable = (DataTable)inputData;
        for(int col=0; col<outputTable.getNColumns(); col++) {
            double[] outputData = outputTable.getColumn(col).getData();
            for(int row=0; row<outputTable.getNRows(); row++) {
                outputData[row] = inputTable.getColumn(row).getData()[col];
            }
        }
        return outputTable;
    }

    /* (non-Javadoc)
     * @see etomica.data.DataProcessor#processDataInfo(etomica.DataInfo)
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        DataTable.Factory factory = (DataTable.Factory)inputDataInfo.getDataFactory();
        int nRows = factory.getNRows();
        int nCols = factory.getNColumns();
        outputTable = new DataTable(inputDataInfo.getLabel(), inputDataInfo.getDimension(), nRows, nCols);
        return outputTable.getDataInfo();
    }

    /* (non-Javadoc)
     * @see etomica.DataSink#getDataCaster(etomica.DataInfo)
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        if(dataInfo.getDataClass() == DataTable.class) {
            return null;
        }
        return new CastToTable();
    }

    private DataTable outputTable;
}
