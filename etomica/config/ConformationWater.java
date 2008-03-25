package etomica.config;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IConformation;
import etomica.space.Space;


public class ConformationWater implements IConformation, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private double bondLengthOH = 4.0;
    private double angleHOH = 109.5*Math.PI/180.;

    public ConformationWater(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomSet list) {
        
        double x = 6.0;
        double y = 6.0;
        
        IAtomPositioned o = (IAtomPositioned)list.getAtom(0);
        o.getPosition().E(new double[] {x, y, 0.0});
               
        IAtomPositioned h1 = (IAtomPositioned)list.getAtom(1);
        h1.getPosition().E(new double[] {x+bondLengthOH, y, 0.0});
                
        IAtomPositioned h2 = (IAtomPositioned)list.getAtom(2);
        h2.getPosition().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});
    }
        
    protected final Space space;
}
