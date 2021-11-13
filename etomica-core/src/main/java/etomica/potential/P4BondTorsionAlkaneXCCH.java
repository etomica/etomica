/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.RandomPositionSourceRectangular;
import etomica.space.*;
import etomica.space3d.Space3D;
import etomica.util.random.IRandom;
import etomica.util.random.RandomNumberGenerator;

/**
 * slight modification from P4BondTorsion class(change energyAtAngle method)
 * Siepmann's Alkane TraPPE-EH model, XCCH torsion potential, X can be H or C, H is for H on CH3 only
 * U(torsion) = Cx * (1 -cos(3phi)); C(C):854K, C(H):717K
 *
 * @author shu
 * Mar 2013
 */
public class P4BondTorsionAlkaneXCCH extends Potential implements PotentialSoft, IPotentialBondTorsion {
    public P4BondTorsionAlkaneXCCH(Space space, double a0, double a1, double a2, double a3) {

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

        gradient = new Vector[4];
        for (int i=0; i<4; i++) {
            gradient[i] = space.makeVector();
        }
    }

    public double energy(IAtomList atomSet) {
        IAtom atom0 = atomSet.get(0);
        IAtom atom1 = atomSet.get(1);
        IAtom atom2 = atomSet.get(2);
        IAtom atom3 = atomSet.get(3);
        dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
        dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
        
        boundary.nearestImage(dr21);
        boundary.nearestImage(dr23);
        boundary.nearestImage(dr34);
        
        double dr23Sq = dr23.squared();
        dr21.PEa1Tv1(-dr21.dot(dr23)/dr23Sq, dr23);
        dr34.PEa1Tv1(-dr34.dot(dr23)/dr23Sq, dr23);

        double cosphi = dr21.dot(dr34) / Math.sqrt(dr21.squared() * dr34.squared());
        // :::::::::::::: check torsion angles ::::::::::::::::: //
//        	int i0 = atom0.getIndex();
//        	int i1 = atom1.getIndex();
//        	int i2 = atom2.getIndex();
//        	int i3 = atom3.getIndex();
//        		System.out.println(String.format("%2d %2d %2d %2d %f", i0, i1, i2, i3, cosphi));
        return u(cosphi);

    }

    @Override
    public double u(double cosphi) {
        double cos2phi = 2 * cosphi * cosphi - 1;
        double cos3phi = cosphi * (2 * cos2phi - 1);
        return a0 + a1 * (1 + cosphi) + a2 * (1 - cos2phi) + a3 * (1 - cos3phi);
        //original:   return a0 + a1*(1+cosphi) + a2*(1-cos2phi) + a3*(1+cos3phi);

    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        double cos2theta = 2 * costheta * costheta - 1;
        double cos3theta = costheta * (2 * cos2theta - 1);
        u[0] = a0 + a1 * (1 + costheta) + a2 * (1 - cos2theta) + a3 * (1 - cos3theta);
        du[0] = 12.0 * a3 * cos2theta - 4.0 * a2 * costheta + a1 - 3 * a3;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public Vector[] gradient(IAtomList atoms) {
        IAtom atom0 = atoms.get(0);
        IAtom atom1 = atoms.get(1);
        IAtom atom2 = atoms.get(2);
        IAtom atom3 = atoms.get(3);
        dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
        dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
        
        boundary.nearestImage(dr21);
        boundary.nearestImage(dr23);
        boundary.nearestImage(dr34);
        
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

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    protected final Vector dr21, dr23, dr34;
    protected final Vector v1, v2;
    protected final Vector gtmp;
    protected Boundary boundary;
    protected double a0, a1, a2, a3;
    protected final Vector[] gradient;
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P4BondTorsionAlkaneXCCH potential = new P4BondTorsionAlkaneXCCH(space, 0, 0, 0, 30);
        IRandom random = new RandomNumberGenerator();
        Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
        RandomPositionSourceRectangular positionSource = new RandomPositionSourceRectangular(space, random);
        positionSource.setBox(box);
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        Atom atom3 = new Atom(space);
        AtomArrayList atoms = new AtomArrayList(4);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
        int n = 40;
        Vector gradient = space.makeVector();
        Vector dr = space.makeVector();
        for (int i=0; i<n; i++) {
            atom0.getPosition().E(positionSource.randomPosition());
            atom1.getPosition().E(positionSource.randomPosition());
            atom2.getPosition().E(positionSource.randomPosition());
            atom3.getPosition().E(positionSource.randomPosition());
            
            double U = potential.energy(atoms);

            int iRand = random.nextInt(4);
            IAtom atom = atoms.get(iRand);
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
