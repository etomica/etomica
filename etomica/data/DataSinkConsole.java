package etomica.data;
import java.io.PrintStream;

import etomica.DataSink;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.units.Dimension;
import etomica.units.Unit;

/**
 * Writes data to console or another print stream.
 */
public class DataSinkConsole implements DataSink, EtomicaElement {

    /**
     * Makes class using System.out as the default output stream.
     */
    public DataSinkConsole() {
        this(System.out);
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
     * Causes the given values to be written to the print stream.
     * Data are written one value per line, all following a header
     * based on the current label and unit.
     */
    public void putData(double[] data) {
        out.println(label + " (" + unit.toString() + ")");
        for(int i=0; i<data.length; i++) {
            out.println(data[i]);
        }
        out.println();
    }

    /**
     * Sets the dimension of the data written to this sink.
     * Part of the DataSink interface.
     */
    public void setDimension(Dimension dimension) {
        this.dimension = dimension;
        if(unit == Unit.UNDEFINED) unit = dimension.defaultIOUnit();
    }
    
    /**
     * @return Returns the label.
     */
    public String getLabel() {
        return label;
    }
    /**
     * @param label The label to set.
     */
    public void setLabel(String label) {
        this.label = label;
    }
    /**
     * @return Returns the unit.
     */
    public Unit getUnit() {
        return unit;
    }
    /**
     * @param unit The unit to set.
     */
    public void setUnit(Unit unit) {
        this.unit = unit;
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
    
    private String label;
    private Dimension dimension = Dimension.UNDEFINED;
    private Unit unit = Unit.UNDEFINED;
    private PrintStream out = System.out;
}