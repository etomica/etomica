package etomica.log;
import etomica.*;
import java.io.*;
import java.util.*;
import java.text.*;

/**
 * Abstract class, used to define the basic methods of reading in/writing out the data.
 * Obtain the status from the integrator(START, INITIALIZE, INTERVAL, DONE),
 * and write() at the specified INTERVAL.
 * 4/15/03
 */

public abstract class LoggerAbstract extends SimulationElement implements Integrator.IntervalListener,
                                                                            java.io.Serializable {
    
    public LoggerAbstract(){
        this(Simulation.instance);
    }
    
    public LoggerAbstract(Simulation sim, Integrator integrator){
        this(sim);
        setIntegrator(integrator);//for a given integrator.
    }
    
    public LoggerAbstract(Simulation sim){
        super(sim, LoggerAbstract.class);
        setUpdateInterval(100);
    }
    
    
    /**
     * Identifies the integrator that fires interval events causing the logger
     * to write to file.  If the integrator is already defined, remove this as
     * listener to it; add this as interval listener to new integrator.  If new
     * integrator is null, removes as listener to current integrator and sets
     * current integrator to null.  Will also add/remove as collision listener
     * if subclass implments IntegratorHard.CollisionListener interface.
     */
    public void setIntegrator(Integrator integratorNew) {
        if(integrator != null){
            integrator.removeIntervalListener(this);
            if(integrator instanceof IntegratorHard && this instanceof IntegratorHard.CollisionListener) 
	            ((IntegratorHard)integrator).removeCollisionListener((IntegratorHard.CollisionListener)this);
	    }
	    
        this.integrator = integratorNew;
        if(integrator == null) return;
        
        integrator.addIntervalListener(this);
        if(integrator instanceof IntegratorHard && this instanceof IntegratorHard.CollisionListener) 
	            ((IntegratorHard)integrator).addCollisionListener((IntegratorHard.CollisionListener)this);
    }
        
    public Integrator getIntegrator() {return integrator;}
    
    
    /**
     * Call write() at the specified INTERVAL.
     */
    public void intervalAction(Integrator.IntervalEvent evt){
        //if(evt.type() == Integrator.IntervalEvent.START) openFile(); //open a file when START.
        if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //only act on INTERVAL events.
	    if(--count == 0) {
	        count = updateInterval;
	        openFile();
	        try {
		        write();
	        } catch(java.io.IOException ex) {
	        	System.out.println("Exception when writing to log ");
	        	ex.printStackTrace();
	        }
	        if(closeFileEachTime) closeFile();
	    }
        if(evt.type() == Integrator.IntervalEvent.DONE) closeFile(); //close the file when DONE.
    }
    
    /**
     * This method will be defined by each subclass to do the log things. e.g.
     * write out data/configuration.
     */
    	//EXAMPLE
		//	protected void write() throws java.io.IOException {
		//		fileWriter.write(Double.toString(meter.average()));
		//	}
    protected abstract void write() throws java.io.IOException;
    
    /**
     * Define how many number of interval events received between writings to
     * file.
     */
    public void setUpdateInterval(int updateInterval){
        this.updateInterval = updateInterval;
        count = updateInterval;
    }
    
    public int getUpdateInterval() {return updateInterval;}
    
    
    private void openFile(){
        try { 
            if(sameFileEachTime && fileIsOpen) return;
            if(fileName == "") fileName = defaultFileName(); //if fileName is not defined yet, use the current date to be the fileName.
            fileWriter = new FileWriter(fileName + fileNameSuffix, appending);
            fileIsOpen = true;
        }catch(IOException e) {
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
        if(fileName.endsWith(fileNameSuffix)) fileName = fileName.substring(0, fileName.length()-4);
        this.fileName = fileName;
    }
    
    public String getFileName() { return (fileName+fileNameSuffix); }
    
    /**
     * Overwrite or attach the data to the file.  Default is to append, not
     * overwrite.
     */
    public void setAppending(boolean appending) {this.appending = appending; }
    
    public boolean isAppending() {return appending;}
    
    /**
     * Write to the same file at each time INTERVAL, and close the file at the
     * end of the simulation. Default is true.
     */
    public void setSameFileEachTime(boolean sameFileEachTime) { this.sameFileEachTime = sameFileEachTime;}
    
    public boolean isSameFileEachTime() {return sameFileEachTime;}
    
    /**
     * Close the file and open a new one at each time INTERVAL. The new file will be named by defaultFileName().
     */
    public void setCloseFileEachTime(boolean closeFileEachTime) { this.closeFileEachTime = closeFileEachTime;}
    
    public boolean isCloseFileEachTime() {return closeFileEachTime;}
        
    private int count; //counts number of events received before calling doAction().
    private int updateInterval; //number of times specified by the user between two actions.
                                //could be set by setUpdateInterval() method.
    private String fileName; //output file name. This name should not have a suffix.
    private boolean appending = true; //whether to overwrite to the existing data 
                                        //or attach the data. Default is not overwriting.
    private Integrator integrator;
    protected FileWriter fileWriter; //this is a handle that subclass can use to write data.
    protected String fileNameSuffix = ".txt"; //Default suffix is text file, subclass can change this.
    private boolean sameFileEachTime = true; //whether to write to the same file at each INTERVAL.
    private boolean closeFileEachTime = false; //whether to close the file and open a new one at each INTERVAL.
    private boolean fileIsOpen = false; //at the beginning, it is false.
    
}