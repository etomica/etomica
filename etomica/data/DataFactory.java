package etomica.data;

import etomica.units.Dimension;


/**
 * Interface for a class that constructs a Data object.  A DataFactory instance is
 * held by DataInfo, and constructs Data instances of the same type and structure
 * as the Data holding the DataInfo.  The type and structure of the Data built by
 * a DataFactory are set at construction and cannot be changed.
 *
 * @author David Kofke
 * 
 * @see Data
 * @see DataInfo
 */

/*
 * History
 * Created on Jul 9, 2005 by kofke
 */
public interface DataFactory {

    /**
     * Constructs a Data instance with a DataInfo having the given label and dimension fields.
     */
    public Data makeData(String label, Dimension dimension);
    
    /**
     * Returns the Class of the Data made by this factory.
     */
    public Class getDataClass();
    
}
