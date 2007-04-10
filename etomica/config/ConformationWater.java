package etomica.config;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.space.Space;


public class ConformationWater extends Conformation {

    private static final long serialVersionUID = 1L;
    private double bondLengthOH = 4.0;
    private double angleHOH = 109.5*Math.PI/180.;
    private final AtomIteratorArrayListSimple moleculeIterator;

    public ConformationWater(Space space) {
        super(space);
        moleculeIterator = new AtomIteratorArrayListSimple();
    }
    
    public void initializePositions(AtomArrayList list) {
        moleculeIterator.setList(list);
        
        double x = 6.0;
        double y = 6.0;
        
        moleculeIterator.reset();
        
        AtomLeaf o = (AtomLeaf)moleculeIterator.nextAtom();
        o.getPosition().E(new double[] {x, y, 0.0});
               
        AtomLeaf h1 = (AtomLeaf)moleculeIterator.nextAtom();
        h1.getPosition().E(new double[] {x+bondLengthOH, y, 0.0});
                
        AtomLeaf h2 = (AtomLeaf)moleculeIterator.nextAtom();
        h2.getPosition().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

    }
        
}
