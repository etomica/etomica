package etomica.models.propane;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * Conformation for 3-point water molecule.
 */
public class ConformationUAPropane extends Conformation {

    private double bondLength = 1.54;
    private double bondAngle = 114.0*Math.PI/180.;

    public ConformationUAPropane(Space space) {
        super(space);
    }
    
    public void initializePositions(AtomArrayList list){
        
        double x = 0.0;
        double y = 0.0;
        
        AtomLeaf UA1 = (AtomLeaf)list.get(0);
        UA1.getPosition().E(new double[] {x, y, 0.0});
               
        AtomLeaf UA2 = (AtomLeaf)list.get(1);
        UA2.getPosition().E(new double[] {x+bondLength, y, 0.0});
                
        AtomLeaf UA3 = (AtomLeaf)list.get(2);
        UA3.getPosition().E(new double[] {x+bondLength+bondLength*Math.cos(Math.PI-bondAngle), y+bondLength*Math.sin(Math.PI-bondAngle), 0.0});
    }
    
    private static final long serialVersionUID = 1L;
}
