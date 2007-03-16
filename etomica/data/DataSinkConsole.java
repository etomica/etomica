package etomica.data;
import java.io.PrintStream;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.util.NameMaker;

/**
 * Writes data to console or another print stream.
 */
public class DataSinkConsole implements DataSink, EtomicaElement, java.io.Serializable {

    /**
     * Makes class using System.out as the default output stream.
     */
    public DataSinkConsole() {
        this(System.out);
        setName(NameMaker.makeName(this.getClass()));
    }
    
    /**
     * Constructor that permits specification of the output printstream.
     * @param outputStream
     */
    public DataSinkConsole(PrintStream outputStream) {
        this.out = outputStream;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pipes data to console");
        return info;
    }
    
    /**
     * Returns null, indicating that any type of Data can be put here without casting.
     */
    public DataProcessor getDataCaster(IDataInfo dataInfo) {
        return null;
    }
    
    /**
     * Causes the given DataInfo to be written to the print stream.
     */
    public void putDataInfo(IDataInfo dataInfo) {
        out.println(dataInfo.toString());
    }
    
    /**
     * Causes the given values to be written to the print stream.
     */
    public void putData(Data data) {
        out.println(data.toString());
    }

    /**
     * Method called to express incredulity.  Short for
     * "No way! Get out of here!"
     * @return Returns the output printstream to which data is written.  
     */
    public PrintStream getOut() {
        return out;
    }
    
    /**
     * Sets the output print stream where data is written.  Default is
     * System.out
     */
    public void setOut(PrintStream out) {
        this.out = out;
    }
    
    /**
     * @return Returns the name.
     */
    public String getName() {
        return name;
    }
    /**
     * @param name The name to set.
     */
    public void setName(String name) {
        this.name = name;
    }

    private static final long serialVersionUID = 1L;
    private PrintStream out = System.out;
    private String name;
}
