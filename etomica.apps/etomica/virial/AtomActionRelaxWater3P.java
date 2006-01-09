/**
 * 
 */
package etomica.virial;

import etomica.action.AtomActionAdapter;
import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.models.water.AtomTreeNodeWater3P;
import etomica.space.Space;
import etomica.space.Vector;

public class AtomActionRelaxWater3P extends AtomActionAdapter {
    public AtomActionRelaxWater3P(Space space) {
        work = space.makeVector();
        cosAngle = Math.cos(109.5/180.0*Math.PI);
        sinAngle = Math.sin(109.5/180.0*Math.PI);
        distance = 1.0;
    }
    
    public void actionPerformed(Atom molecule) {
        AtomTreeNodeWater3P waterNode = (AtomTreeNodeWater3P)molecule.node;
        AtomLeaf O = waterNode.O;
        AtomLeaf H1 = waterNode.H1;
        AtomLeaf H2 = waterNode.H2;
        // normalize OH1
        Vector p1 = H1.coord.position();
        p1.ME(O.coord.position());
        p1.TE(1/Math.sqrt(p1.squared()));
        Vector p2 = H2.coord.position();
        p2.ME(O.coord.position());
        p2.TE(1/Math.sqrt(p2.squared()));
        // move H2 to fix bond angle
        double d = p1.dot(p2);
        work.Ev1Pa1Tv2(p2,-d,p1);
        work.TE(1/Math.sqrt(work.squared()));
        p2.Ea1Tv1(sinAngle,work);
        p2.PEa1Tv1(cosAngle,p1);
        p2.TE(distance/Math.sqrt(p2.squared()));
        p1.TE(distance);
        p1.PE(O.coord.position());
        p2.PE(O.coord.position());
    }

    private final Vector work;
    private final double sinAngle, cosAngle, distance;
}