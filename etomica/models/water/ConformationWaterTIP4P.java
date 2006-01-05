package etomica.models.water;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * Conformation for 4-point water molecule.
 */
public class ConformationWaterTIP4P extends Conformation {

    private double bondLengthOH = 0.9572;
    private double angleHOH = 104.52*Math.PI/180.;
    private double rOM=0.15;
    private final AtomIteratorArrayListSimple iterator;

    public ConformationWaterTIP4P(Space space) {
        super(space);
        iterator = new AtomIteratorArrayListSimple();
    }
    
    public void initializePositions(AtomArrayList list){
        
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
        
        Atom m = iterator.nextAtom();
        m.coord.position().E(new double[] {x+rOM*Math.cos(angleHOH/2.0), y+rOM*Math.sin(angleHOH/2.0), 0.0});

    }//end of initializePositions
    
    
}
