package etomica.potential;

import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space3d.IVector3D;

/**
 * Torsion potential.
 *    Jorgensen, W. L., Madura, J. D. and Swenson, C. J.
 *    J. Am. chem. Soc., 106, 6638-6646 (1984)
 *
 * @author Andrew Schultz
 */
public class P4BondTorsion extends Potential {

    public P4BondTorsion(Space space, double a1, double a2, double a3) {
        super(4, space);
        dr21 = space.makeVector();
        dr23 = space.makeVector();
        dr34 = space.makeVector();
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;
    }

    public void setBox(Box box) {
        nearestImageTransformer = box.getBoundary();
    }

    public double energy(AtomSet atomSet) {
        IAtomPositioned atom0 = (IAtomPositioned)atomSet.getAtom(0);
        IAtomPositioned atom1 = (IAtomPositioned)atomSet.getAtom(1);
        IAtomPositioned atom2 = (IAtomPositioned)atomSet.getAtom(2);
        IAtomPositioned atom3 = (IAtomPositioned)atomSet.getAtom(3);
        dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
        dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
        
        nearestImageTransformer.nearestImage(dr21);
        nearestImageTransformer.nearestImage(dr23);
        nearestImageTransformer.nearestImage(dr34);
        
        double dr23Sq = dr23.squared();
        dr21.PEa1Tv1(-dr21.dot(dr23)/dr23Sq, dr23);
        dr34.PEa1Tv1(-dr34.dot(dr23)/dr23Sq, dr23);
        
        double cosphi = dr21.dot(dr34)/Math.sqrt(dr21.squared()*dr34.squared());
        return energyAtAngle(cosphi);
    }
    
    public double energyAtAngle(double cosphi) {
        double cos2phi = 2*cosphi*cosphi-1;
        double cos3phi = cosphi*(2*cos2phi-1);

        return a1*(1+cosphi) + a2*(1-cos2phi) + a3*(1+cos3phi);
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    protected final IVector dr21, dr23, dr34;
    protected NearestImageTransformer nearestImageTransformer;
    protected double a1, a2, a3;
    private static final long serialVersionUID = 1L;
}
