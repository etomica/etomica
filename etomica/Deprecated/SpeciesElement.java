package simulate;
import java.util.*; // for neighborList

public interface SpeciesElement
{
    void setRm(double rm);
    
    double getCollisionTime();
    void setCollisionTime(double collisionTime);
    void decrementCollisionTime(double interval);
    
    SpeciesElement getCollisionPartner();
    void setCollisionPartner(SpeciesElement partner);
    
    SpeciesElement getNext();
    void setNext(SpeciesElement element);

    SpeciesElement getPrevious();
    void setPrevious(SpeciesElement element);

    int getSpeciesIndex();
    void setSpeciesIndex(int index);
    
    void zeroForce();
    void addForce(double[] force);
    void subtractForce(double[] force);
    
    double getKineticEnergy();

    public void addNeighbor(SpeciesElement e);
    public void clearNeighborList();
    public Enumeration getNeighborList();
    
}