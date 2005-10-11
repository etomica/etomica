package etomica.data;

public interface DataPipeForked extends DataPipe {

    public int getDataSinkCount();

    /**
     * @return the i-th DataSink
     */
    public DataSink getDataSink(int i);

    /**
     * Sets the list of DataSinks that receive the Data entering this DataFork.
     * All previously added DataSinks are discarded.  If argument is null, all
     * existing DataSinks are discarded and none are added.
     * 
     * @param dataSinks The data sinks to set.
     */
    public void setDataSinks(DataSink[] dataSinks);

    /**
     * Adds the given DataSink to those receiving the Data entering this DataFork,
     * keeping all previously entered DataSinks.  If argument is null, no action
     * is performed.
     * 
     * @param dataSink
     */
    public void addDataSink(DataSink newDataSink);

    /**
     * Removes the specified data sink.  Does nothing if the given DataSink is
     * not currently a DataSink for this DataFork.
     * 
     * @param dataSink data sink to be removed from this list, if present.
     */
    public void removeDataSink(DataSink dataSink);

}