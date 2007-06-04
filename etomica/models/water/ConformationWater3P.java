package etomica.models.water;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtomPositioned;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * Conformation for 3-point water molecule.
 */
public class ConformationWater3P extends Conformation {

    private double bondLengthOH = 1.0;
    private double angleHOH = 109.5*Math.PI/180.;

    public ConformationWater3P(Space space) {
        super(space);
    }
    
    public void initializePositions(AtomArrayList list){
        
        double x = 0.0;
        double y = 0.0;
        
        IAtomPositioned o = (IAtomPositioned)list.getAtom(0);
        o.getPosition().E(new double[] {x, y, 0.0});
               
        IAtomPositioned h1 = (IAtomPositioned)list.getAtom(0);
        h1.getPosition().E(new double[] {x+bondLengthOH, y, 0.0});
                
        IAtomPositioned h2 = (IAtomPositioned)list.getAtom(0);
        h2.getPosition().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});

    }
    
    private static final long serialVersionUID = 1L;
}
