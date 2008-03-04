package etomica.models.water;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * Conformation for 3-point water molecule.
 */
public class ConformationWater3P extends Conformation {

    protected static final double bondLengthOH = 1.0;
    protected static final double angleHOH = 109.5*Math.PI/180.;

    public ConformationWater3P(Space space) {
        super(space);
    }
    
    public void initializePositions(IAtomSet list){
        
        IAtomPositioned o = (IAtomPositioned)list.getAtom(2);
        o.getPosition().E(new double[] {0, 0, 0.0});

        double x = bondLengthOH*Math.sin(0.5*angleHOH);
        double y = bondLengthOH*Math.cos(0.5*angleHOH);
        
        IAtomPositioned h1 = (IAtomPositioned)list.getAtom(0);
        h1.getPosition().E(new double[] {-x, y, 0.0});
                
        IAtomPositioned h2 = (IAtomPositioned)list.getAtom(1);
        h2.getPosition().E(new double[] {+x, y, 0.0});
    }
    
    private static final long serialVersionUID = 1L;
}
