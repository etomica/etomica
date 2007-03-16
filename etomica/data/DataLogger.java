package etomica.data;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;

/**
 * DataSink that manages a FileWriter and also listens to non-interval 
 * integrator events and sends appropriate data to a DataWriter.
 */
public class DataLogger implements DataPipe, IntegratorNonintervalListener, java.io.Serializable {
    
    public DataLogger(){
        setWriteInterval(100);
        setPriority(300);
        tag = new DataTag();
    }
    
    /**
     * Gives data to DataSink for writing
     */
    public void putData(Data data) {
        if (--count == 0) {
            count = writeInterval;
            doWrite(data);
        }
    }
    
    public void putDataInfo(IDataInfo newDataInfo) {
        dataInfo = newDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        insertTransformerIfNeeded();
        if (dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    protected void insertTransformerIfNeeded() {
        if(dataWriter == null || dataInfo == null) return;
        //remove transformer if one was previously inserted
        dataSink = dataWriter;
        DataProcessor caster = dataWriter.getDataCaster(dataInfo);
        if(caster != null) {
            caster.setDataSink(dataWriter);
            dataSink = caster;
        }
        dataSink.putDataInfo(dataInfo);
    }

    public DataProcessor getDataCaster(IDataInfo dataInfo) {
        // we don't care about the type although the DataWriter might
        return null;
    }
    
    /**
     * Sets the data writer to the given DataSink.  It must be an instance of
     * DataLogger.DataWriter
     */
    public void setDataSink(DataSink dataWriter) {
        this.dataWriter = (DataWriter)dataWriter;
        insertTransformerIfNeeded();
    }
    
    public DataSink getDataSink() {
        return dataWriter;
    }
    
    /**
     * Close file when integrator is done.
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt){
        if(evt.type() == IntegratorNonintervalEvent.DONE && dataWriter != null) {
            if (writeOnFinishSource != null) {
                doWrite(writeOnFinishSource.getData());
            }
            closeFile(); //close the file when DONE.
        }
    }
    
    /**
     * Performs the open/write/close file actions in accordance with the
     * settings of the Logger.  Called by intervalAction and by the
     * actionPerformed method of any Actions made by the writeAction method.
     */
    private void doWrite(Data data) {
        openFile();
        dataWriter.setFileWriter(fileWriter);
        dataSink.putData(data);
        if(closeFileEachTime || !appending) {
            closeFile();
        }
    }
    
    /**
     * Define how many number of interval events received between writings to
     * file.
     */
    public void setWriteInterval(int newWriteInterval){
        writeInterval = newWriteInterval;
        count = writeInterval;
    }
    
    public int getWriteInterval() {
        return writeInterval;
    }
    
    
    private void openFile(){
        try { 
            if(sameFileEachTime && fileIsOpen) return;
            if(fileName == "") fileName = defaultFileName(); //if fileName is not defined yet, use the current date to be the fileName.
            fileWriter = new FileWriter(fileName + fileNameSuffix, appending);
            if (!appending) {
                dataWriter.reset();
            }
            fileIsOpen = true;
        }
        catch(IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
    }
    
    protected String defaultFileName() {//now use the current date to be the fileName, subclass could change this.
        SimpleDateFormat formatter = new SimpleDateFormat("MMMd_H_mm_ss", new Locale("en","US"));
        Date now = new Date();
        String currentTime = formatter.format(now);
        return currentTime;
    }
    
    private void closeFile(){ 
        try {
            if(fileIsOpen) { fileWriter.close(); fileIsOpen = false;}
            if(!sameFileEachTime) fileName = "";
        }catch(IOException e) {
            System.err.println("Cannot close a file, caught IOException: " + e.getMessage());
        }
    }
    
    /**
     * Define the name of the output file.
     */
    public void setFileName(String fileName) {//test if the fileName ends with ".txt", cut ".txt" from it.
        if(fileName.endsWith(fileNameSuffix)) {
            fileName = fileName.substring(0, fileName.length()-4);
        }
        this.fileName = fileName;
    }
    
    public String getFileName() {
        return (fileName+fileNameSuffix);
    }
    
    /**
     * Overwrite or attach the data to the file.  Default is to append, not
     * overwrite.
     */
    public void setAppending(boolean appending) {
        this.appending = appending; 
    }
    
    public boolean isAppending() {
        return appending;
    }
    
    /**
     * Write to the same file at each time INTERVAL, and close the file at the
     * end of the simulation. Default is true.
     */
    public void setSameFileEachTime(boolean sameFileEachTime) {
        this.sameFileEachTime = sameFileEachTime;
    }
    
    public boolean isSameFileEachTime() {
        return sameFileEachTime;
    }
    
    /**
     * Close the file and open a new one at each time INTERVAL. The new file will be named by defaultFileName().
     */
    public void setCloseFileEachTime(boolean closeFileEachTime) {
        this.closeFileEachTime = closeFileEachTime;
    }
    
    public boolean isCloseFileEachTime() {
        return closeFileEachTime;
    }
    
    /**
     * @return Returns the writeOnFinish.
     */
    public DataSource getWriteOnFinishDataSource() {
        return writeOnFinishSource;
    }

    /**
     * @param writeOnFinish The writeOnFinish to set.
     */
    public void setWriteOnFinish(DataSource newWriteOnFinishSource) {
        writeOnFinishSource = newWriteOnFinishSource;
    }

    /**
     * @return Returns the writeOnInterval.
     */
    public boolean isWriteOnInterval() {
        return writeOnInterval;
    }

    /**
     * @param writeOnInterval The writeOnInterval to set.
     */
    public void setWriteOnInterval(boolean writeOnInterval) {
        this.writeOnInterval = writeOnInterval;
    }

    /**
     * @return Returns the interval-listener priority.
     */
    public int getPriority() {
        return priority;
    }
    
    /**
     * Sets the interval-listener priority.  Default value is 300, which
     * puts this after central-image enforcement and accumulator updates.
     * @param priority The priority to set.
     */
    public void setPriority(int priority) {
        this.priority = priority;
    }


    private static final long serialVersionUID = 2L;
    private int count; //counts number of events received before calling doAction().
    private int writeInterval; //number of times specified by the user between two actions.
                                //could be set by setUpdateInterval() method.
    private String fileName; //output file name. This name should not have a suffix.
    private boolean appending = true; //whether to overwrite to the existing data 
                                        //or attach the data. Default is not overwriting.
    protected transient FileWriter fileWriter; //this is a handle that subclass can use to write data.
    protected String fileNameSuffix = ".txt"; //Default suffix is text file, subclass can change this.
    private boolean sameFileEachTime = true; //whether to write to the same file at each INTERVAL.
    private boolean closeFileEachTime = false; //whether to close the file and open a new one at each INTERVAL.
    private transient boolean fileIsOpen = false; //at the beginning, it is false.
    private int priority;
    private DataSource writeOnFinishSource = null;
    private boolean writeOnInterval = true;
    private DataWriter dataWriter;
    private DataSink dataSink;
    private IDataInfo dataInfo;
    protected final DataTag tag;

    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        in.defaultReadObject();
        fileWriter = null;
        // the file isn't open!  :)
        // this will also cause fileWriter to be recreated when needed
        fileIsOpen = false;
    }
    
    /**
     * Interface for a DataSink that actually writes data to a file
     */
    public interface DataWriter extends DataSink {
        
        /**
         * Sets the FileWriter to be used for actual file I/O.
         */
        public void setFileWriter(FileWriter newFileWriter);

        /**
         * Informs the DataSink that it will be writing the beginning of the 
         * file so printing header information is appropriate.
         */
        public void reset();
    }
}
