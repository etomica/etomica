/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.action.activity.ControllerEvent;
import etomica.util.IEvent;
import etomica.util.IListener;

import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

/**
 * DataSink that manages a FileWriter and also listens to non-interval 
 * integrator events and sends appropriate data to a DataWriter.
 */
public class DataLogger extends DataProcessor implements IListener, java.io.Serializable {

    private static final long serialVersionUID = 2L;
    //or attach the data. Default is not overwriting.
    protected transient FileWriter fileWriter; //this is a handle that subclass can use to write data.
    protected String fileNameSuffix = ".dat"; //Default suffix is text file, subclass can change this.
    private int count; //counts number of events received before calling doAction().
    private int writeInterval; //number of times specified by the user between two actions.
    //could be set by setUpdateInterval() method.
    private String fileName; //output file name. This name should not have a suffix.
    private boolean appending = true; //whether to overwrite to the existing data
    private boolean sameFileEachTime = true; //whether to write to the same file at each INTERVAL.
    private boolean closeFileEachTime = false; //whether to close the file and open a new one at each INTERVAL.
    private transient boolean fileIsOpen = false; //at the beginning, it is false.
    private IEtomicaDataSource writeOnFinishSource = null;
    private boolean writeOnInterval = true;
    
    public DataLogger(){
        super();
        setWriteInterval(100);
    }
    
    public void setDataSink(IDataSink dataSink) {
        if (!(dataSink instanceof DataWriter)) {
            throw new IllegalArgumentException("data sink must be a DataWriter");
        }
        super.setDataSink(dataSink);
    }
    
    /**
     * Gives data to DataSink for writing
     */
    public IData processData(IData data) {

        return data;
    }
    
    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo newDataInfo) {
        dataInfo = newDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        return dataInfo;
    }
    
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        // we don't care about the type although the DataWriter might
        return null;
    }
    
    /**
     * Close file when integrator is done.
     */
    public void actionPerformed(IEvent evt){
        if(evt instanceof ControllerEvent) {
            if(dataSink != null && (((ControllerEvent)evt).getType() == ControllerEvent.NO_MORE_ACTIONS ||
                    ((ControllerEvent)evt).getType() == ControllerEvent.HALTED)) {
                if (writeOnFinishSource != null) {
                    putData(writeOnFinishSource.getData());
                }
                closeFile(); //close the file when finished
            }
        }
    }
    
    /**
     * Performs the open/write/close file actions in accordance with the
     * settings of the Logger.  Called by intervalAction and by the
     * actionPerformed method of any Actions made by the writeAction method.
     */
    public void putData(IData data) {
        openFile();
        ((DataWriter)trueDataSink).setFileWriter(fileWriter);
        super.putData(data);
        if (closeFileEachTime) {
            closeFile();
        }
    }

    public int getWriteInterval() {
        return writeInterval;
    }
    
    /**
     * Define how many number of interval events received between writings to
     * file.
     */
    public void setWriteInterval(int newWriteInterval){
        writeInterval = newWriteInterval;
        count = writeInterval;
    }

    private void openFile(){
        try {
            if(sameFileEachTime && fileIsOpen) return;
            if(fileName == "") fileName = defaultFileName(); //if fileName is not defined yet, use the current date to be the fileName.
            fileWriter = new FileWriter(fileName + fileNameSuffix, appending);
            if (!appending) {
                ((DataWriter)trueDataSink).reset();
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

    public void closeFile() {
        try {
            if(fileIsOpen) { fileWriter.close(); fileIsOpen = false;}
            if(!sameFileEachTime) fileName = "";
        }catch(IOException e) {
            System.err.println("Cannot close a file, caught IOException: " + e.getMessage());
        }
    }

    public String getFileName() {
        return (fileName + fileNameSuffix);
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

    public boolean isAppending() {
        return appending;
    }

    /**
     * Overwrite or attach the data to the file.  Default is to append, not
     * overwrite.
     */
    public void setAppending(boolean appending) {
        this.appending = appending;
    }

    public boolean isSameFileEachTime() {
        return sameFileEachTime;
    }

    /**
     * Write to the same file at each time INTERVAL, and close the file at the
     * end of the simulation. Default is true.
     */
    public void setSameFileEachTime(boolean sameFileEachTime) {
        this.sameFileEachTime = sameFileEachTime;
    }

    public boolean isCloseFileEachTime() {
        return closeFileEachTime;
    }

    /**
     * Close the file and open a new one at each time INTERVAL. The new file will be named by defaultFileName().
     */
    public void setCloseFileEachTime(boolean closeFileEachTime) {
        this.closeFileEachTime = closeFileEachTime;
    }

    /**
     * @return Returns the writeOnFinish.
     */
    public IEtomicaDataSource getWriteOnFinishDataSource() {
        return writeOnFinishSource;
    }

    /**
     * @param newWriteOnFinishSource The writeOnFinish to set.
     */
    public void setWriteOnFinish(IEtomicaDataSource newWriteOnFinishSource) {
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
    public interface DataWriter extends IDataSink {
        
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
