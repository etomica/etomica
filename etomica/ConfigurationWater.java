package etomica;

import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;


public class ConfigurationWater extends Conformation {

    private double bondLengthOH = 4.0;
    private double angleHOH = 109.5*Math.PI/180.;
    private final AtomIteratorListSimple moleculeIterator;

    public ConfigurationWater(Space space) {
        super(space);
        moleculeIterator = new AtomIteratorListSimple();
    }
    
    public void initializePositions(AtomList list) {
        moleculeIterator.setList(list);
        
        double x = 6.0;
        double y = 6.0;
        
        moleculeIterator.reset();
        
        Atom o = moleculeIterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0.0});
               
        Atom h1 = moleculeIterator.nextAtom();
        h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
        Atom h2 = moleculeIterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

    }
        
}
