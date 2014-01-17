package etomica.atom;

import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.units.BohrRadius;

public class AtomHydrogen extends AtomOriented {    
    protected double bondLength;

    public AtomHydrogen(ISpace space, IAtomTypeOriented atype, double bl) {
        super(space, atype);        
//        bondLength = BohrRadius.UNIT.toSim(1.401065676);
        bondLength = bl;//BohrRadius.UNIT.toSim(1.448736);
    }
    public double getBondLength() {
        return bondLength;        
    }
    public void setBondLength(double x) {
        bondLength = x;
    }

}
