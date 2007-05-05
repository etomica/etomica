/**
 * 
 */
package etomica.virial;

import etomica.action.AtomActionAdapter;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.models.water.AtomWater3P;
import etomica.space.IVector;
import etomica.space.Space;

public class AtomActionRelaxWater3P extends AtomActionAdapter {
    public AtomActionRelaxWater3P(Space space) {
        work = space.makeVector();
        cosAngle = Math.cos(109.5/180.0*Math.PI);
        sinAngle = Math.sin(109.5/180.0*Math.PI);
        distance = 1.0;
    }
    
    public void actionPerformed(IAtom molecule) {
        IAtomPositioned O = ((AtomWater3P)molecule).O;
        IAtomPositioned H1 = ((AtomWater3P)molecule).H1;
        IAtomPositioned H2 = ((AtomWater3P)molecule).H2;
        // normalize OH1
        IVector p1 = H1.getPosition();
        p1.ME(O.getPosition());
        p1.TE(1/Math.sqrt(p1.squared()));
        IVector p2 = H2.getPosition();
        p2.ME(O.getPosition());
        p2.TE(1/Math.sqrt(p2.squared()));
        // move H2 to fix bond angle
        double d = p1.dot(p2);
        work.E(p2);
        work.PEa1Tv1(-d,p1);
        work.TE(1/Math.sqrt(work.squared()));
        p2.Ea1Tv1(sinAngle,work);
        p2.PEa1Tv1(cosAngle,p1);
        p2.TE(distance/Math.sqrt(p2.squared()));
        p1.TE(distance);
        p1.PE(O.getPosition());
        p2.PE(O.getPosition());
    }

    private static final long serialVersionUID = 1L;
    private final IVector work;
    private final double sinAngle, cosAngle, distance;
}