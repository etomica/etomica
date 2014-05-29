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
    
    public static double getAvgBondLength(double x) {
        double a0 = 0.0766111;
        double a1 = 6.84704E-07;
        double a2 = -3.89889E-11;
        double a3 = 6.73971E-14;
        double y = 10*(a0 + a1*x + a2*x*x + a3*x*x*x);
//        System.out.println(x+" "+y);
//        System.exit(1);
        return y;
    }

}
