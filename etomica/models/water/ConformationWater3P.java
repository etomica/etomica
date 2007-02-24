package etomica.models.water;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * Conformation for 3-point water molecule.
 */
public class ConformationWater3P extends Conformation {

    private double bondLengthOH = 1.0;
    private double angleHOH = 109.5*Math.PI/180.;
    private final AtomIteratorArrayListSimple iterator;

    public ConformationWater3P(Space space) {
        super(space);
        iterator = new AtomIteratorArrayListSimple();
    }
    
    public void initializePositions(AtomArrayList list){
        
        iterator.setList(list);
        double x = 0.0;
        double y = 0.0;
        
        iterator.reset();
        
        AtomLeaf o = (AtomLeaf)iterator.nextAtom();
        o.getCoord().getPosition().E(new double[] {x, y, 0.0});
               
        AtomLeaf h1 = (AtomLeaf)iterator.nextAtom();
        h1.getCoord().getPosition().E(new double[] {x+bondLengthOH, y, 0.0});
                
        AtomLeaf h2 = (AtomLeaf)iterator.nextAtom();
        h2.getCoord().getPosition().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

    }
    
    private static final long serialVersionUID = 1L;
}
