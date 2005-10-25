package etomica.space;

public class CoordinateFactoryNull implements CoordinateFactory, java.io.Serializable {

    public ICoordinate makeCoordinate() {
        return null;
    }
    
    public void setKinetic(boolean kinetic) {}
    
    public boolean isKinetic() {
        return false;
    }
}
