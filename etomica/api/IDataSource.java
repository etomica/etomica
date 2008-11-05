package etomica.api;

public interface IDataSource {

    /**
     * @return the data given by this source
     */
    public IData getData();

    /**
     * Returns the DataInfo instance that will be held by Data
     * given by this source.  This information is useful for
     * setting up the data stream and for providing annotation when
     * displaying or writing the Data.
     */
//    public IDataInfo getDataInfo();

}
