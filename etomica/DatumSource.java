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
    
    public interface MultiUser {
        public void setDatumSources(DatumSource[] source);
        public void setDatumSources(int i, DatumSource source);
        public DatumSource[] getDatumSources();
        public DatumSource getDatumSources(int i);
        public void addDatumSources(DatumSource source);
        public void addDatumSources(DatumSource[] source);
    }
    
}