package etomica.models.water;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * Conformation for 3-point water molecule.
 */
public class ConformationWater3P extends Conformation {

    private double bondLengthOH = 1.0;
    private double angleHOH = 109.5*Math.PI/180.;
    private final AtomIteratorListSimple iterator;

    public ConformationWater3P(Space space) {
        super(space);
        iterator = new AtomIteratorListSimple();
    }
    
    public void initializePositions(AtomList list){
        
        iterator.setList(list);
        double x = 0.0;
        double y = 0.0;
        
        iterator.reset();
        
        Atom o = iterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0.0});
               
        Atom h1 = iterator.nextAtom();
        h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
        Atom h2 = iterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

    }//end of initializePositions
    
    
}
