package etomica;

/**
 * Interface for a class that can provide a single datum of type double.
 */
public interface DatumSource {
    
    public double value(DataSource.ValueType type);
    
    public String getLabel();
    
    public etomica.units.Dimension getDimension();
        
    public interface User {
        public void setDatumSource(DatumSource source);
        public DatumSource getDatumSource();
    }
}