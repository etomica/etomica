package etomica.space;

/**
 * Returns null for a Coordinate, which is appropriate for Atom groups, which 
 * holds a null coord field.
 */
public class CoordinateFactoryNull implements CoordinateFactory, java.io.Serializable {

    public ICoordinate makeCoordinate() {
        return null;
    }
    
    public void setKinetic(boolean kinetic) {}
    
    public boolean isKinetic() {
        return false;
    }
}
