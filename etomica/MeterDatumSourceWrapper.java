package etomica;

/**
 * Meter that contains another DatumSource object (typically other than a Meter)
 * and which obtains its "measured" values from it.  Useful to apply history
 * or histogram functionality to the DatumSource (such as a Modulator).
 *
 * @author David Kofke
 */
 
 /* See PistonCylinder module for an example, where user movement of pressure slider
  * is recorded this way.
  */
  
 public class MeterDatumSourceWrapper extends Meter {
    
    private DatumSource source;
    private DataSource.ValueType whichValue = null;
    
    public MeterDatumSourceWrapper(DatumSource source) {
        this(Simulation.instance, source);
    }
    
    public MeterDatumSourceWrapper(Simulation sim, DatumSource source) {
        super(sim);
        this.source = source;
        setLabel(source.getLabel());
    }
    
    public DatumSource getDatumSource() {return source;}
    
    /**
     * Accessor method of which value to obtain from the internal DatumSource.
     * Default is null.
     */
    public DataSource.ValueType getWhichValue() {return whichValue;}
    /**
     * Mutator method of which value to obtain from the internal DatumSource.
     * Default is null.
     */
    public void setWhichValue(DataSource.ValueType type) {whichValue = type;}
    
    /**
     * Dimension as prescribed by the internal DatumSource.
     */
    public etomica.units.Dimension getDimension() {return source.getDimension();}
    
    /**
     * Returns the output of a call to the value method of the internal DatumSource,
     * given the whichValue type as an argument.
     */
    public double currentValue() {
        return (source != null) ? source.value(whichValue) : Double.NaN;
    }
    
 }//end of MeterDatumSourceWrapper