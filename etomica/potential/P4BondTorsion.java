package etomica.potential;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.INearestImageTransformer;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.box.Box;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.IVectorRandom;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.util.RandomNumberGenerator;

/**
 * Torsion potential.
 *    Jorgensen, W. L., Madura, J. D. and Swenson, C. J.
 *    J. Am. chem. Soc., 106, 6638-6646 (1984)
 *
 * @author Andrew Schultz
 */
public class P4BondTorsion extends Potential implements PotentialSoft {

    public P4BondTorsion(Space space, double a0, double a1, double a2, double a3) {
        super(4, space);
        dr21 = space.makeVector();
        dr23 = space.makeVector();
        dr34 = space.makeVector();
        this.a0 = a0;
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;
        v1 = space.makeVector();
        v2 = space.makeVector();

        gtmp = space.makeVector();

        gradient = new IVector[4];
        for (int i=0; i<4; i++) {
            gradient[i] = space.makeVector();
        }
    }

    public void setBox(IBox box) {
        nearestImageTransformer = box.getBoundary();
    }

    public double energy(IAtomSet atomSet) {
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

        return a0 + a1*(1+cosphi) + a2*(1-cos2phi) + a3*(1+cos3phi);
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public IVector[] gradient(IAtomSet atoms) {
        IAtomPositioned atom0 = (IAtomPositioned)atoms.getAtom(0);
        IAtomPositioned atom1 = (IAtomPositioned)atoms.getAtom(1);
        IAtomPositioned atom2 = (IAtomPositioned)atoms.getAtom(2);
        IAtomPositioned atom3 = (IAtomPositioned)atoms.getAtom(3);
        dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
        dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
        
        nearestImageTransformer.nearestImage(dr21);
        nearestImageTransformer.nearestImage(dr23);
        nearestImageTransformer.nearestImage(dr34);
        
        double dr23Sq = dr23.squared();
        double dr23dotdr21odr23Sq = dr23.dot(dr21)/dr23Sq;
        double dr23dotdr34odr23Sq = dr23.dot(dr34)/dr23Sq;
        
        v1.E(dr21);
        v1.PEa1Tv1(-dr21.dot(dr23)/dr23Sq, dr23);
        v2.E(dr34);
        v2.PEa1Tv1(-dr23dotdr34odr23Sq, dr23);

        double v1Sq = v1.squared();
        double v2Sq = v2.squared();
        double v1dotv2 = v1.dot(v2);
        double v1v2 = Math.sqrt(v1Sq*v2Sq);
        if (v1v2/dr23Sq < 1e-7) {
            // either 123 or 234 are nearly colinear; evaluating gradient is
            // hard.  let's go shopping
            gradient[0].E(0);
            gradient[1].E(0);
            gradient[2].E(0);
            gradient[3].E(0);
            return gradient;
        }
        double v1v2_3 = v1v2*v1v2*v1v2;
        
        double cosphi = v1.dot(v2)/Math.sqrt(v1Sq*v2Sq);
        double cos2phi = cosphi*cosphi;  // note, this is different than cos2phi in energy()
        
        double dUdcosphi = 12.0*a3*cos2phi - 4.0*a2*cosphi + a1 - 3*a3;

        gradient[0].Ea1Tv1(1.0/v1v2, v2);
        gradient[0].PEa1Tv1(-v1dotv2*v2Sq/v1v2_3, v1);
        gradient[0].TE(dUdcosphi);
        
        // d(v1dotv2)/dr1
        gtmp.Ev1Pv2(dr21, dr23);
        gtmp.TE(dr23dotdr34odr23Sq);
        gtmp.ME(dr34);
        gtmp.PEa1Tv1(dr21.dot(dr23)/dr23Sq, dr34);
        gtmp.PEa1Tv1(-2*dr21.dot(dr23)*dr23dotdr34odr23Sq/dr23Sq, dr23);
        gtmp.TE(1.0/v1v2);
        gradient[1].E(gtmp);

        // d(v1^2)/dr1
        gtmp.Ev1Pv2(dr21, dr23);
        gtmp.TE(2.0*dr21.dot(dr23)/dr23Sq);
        gtmp.PEa1Tv1(-2.0, dr21);
        gtmp.PEa1Tv1(-2.0*dr23dotdr21odr23Sq*dr23dotdr21odr23Sq, dr23);
        gtmp.TE(-0.5*v1dotv2*v2Sq/v1v2_3);
        gradient[1].PE(gtmp);
        
        // d(v2^2)/dr1
        gtmp.Ea1Tv1(2*dr23dotdr34odr23Sq, dr34);
        gtmp.PEa1Tv1(-2*dr23dotdr34odr23Sq*dr23dotdr34odr23Sq, dr23);
        gtmp.TE(-0.5*v1dotv2*v1Sq/v1v2_3);
        gradient[1].PE(gtmp);
        gradient[1].TE(dUdcosphi);

        gradient[3].Ea1Tv1(1.0/v1v2, v1);
        gradient[3].PEa1Tv1(-v1dotv2*v1Sq/v1v2_3, v2);
        gradient[3].TE(dUdcosphi);
        
        gradient[2].Ea1Tv1(-1, gradient[0]);
        gradient[2].PEa1Tv1(-1, gradient[1]);
        gradient[2].PEa1Tv1(-1, gradient[3]);

        return gradient;
    }

    public IVector[] gradient(IAtomSet atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomSet atoms) {
        return 0;
    }

    private static final long serialVersionUID = 1L;
    protected final IVector dr21, dr23, dr34;
    protected final IVector v1, v2;
    protected final IVector gtmp;
    protected INearestImageTransformer nearestImageTransformer;
    protected double a0, a1, a2, a3;
    protected final IVector[] gradient;
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P4BondTorsion potential = new P4BondTorsion(space, 0, 10, 20, 30);
        IRandom random = new RandomNumberGenerator();
        Box box = new Box(new BoundaryRectangularNonperiodic(space, random), space);
        potential.setBox(box);
        AtomLeaf atom0 = new AtomLeaf(space);
        AtomLeaf atom1 = new AtomLeaf(space);
        AtomLeaf atom2 = new AtomLeaf(space);
        AtomLeaf atom3 = new AtomLeaf(space);
        AtomArrayList atoms = new AtomArrayList(4);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
        int n = 40;
        IVector gradient = space.makeVector();
        IVectorRandom dr = (IVectorRandom)space.makeVector();
        for (int i=0; i<n; i++) {
            atom0.getPosition().E(box.getBoundary().randomPosition());
            atom1.getPosition().E(box.getBoundary().randomPosition());
            atom2.getPosition().E(box.getBoundary().randomPosition());
            atom3.getPosition().E(box.getBoundary().randomPosition());
            
            double U = potential.energy(atoms);

            int iRand = random.nextInt(4);
            IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(iRand);
            gradient.E(potential.gradient(atoms)[iRand]);
            
            dr.setRandomSphere(random);
            dr.TE(0.0001);
            double expectedDeltaU = gradient.dot(dr);
            
            atom.getPosition().PE(dr);
            
            double newU = potential.energy(atoms);
            
            System.out.println(expectedDeltaU+" "+(newU-U)+" "+(expectedDeltaU-newU+U));
        }
    }
}
