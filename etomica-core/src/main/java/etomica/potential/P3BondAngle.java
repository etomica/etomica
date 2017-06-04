/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.box.RandomPositionSourceRectangular;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.units.Angle;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.util.RandomNumberGenerator;

/**
 * Simple 3-body soft bond-angle potential 
 * @author andrew
 */
public class P3BondAngle extends Potential implements PotentialSoft {

    public P3BondAngle(Space space) {
        super(3, space);
        dr12 = space.makeVector();
        dr23 = space.makeVector();
        setAngle(Math.PI);
        gradient = new Vector[3];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        gradient[2] = space.makeVector();
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IAtomList atomSet) {
        IAtom atom0 = atomSet.getAtom(0);
        IAtom atom1 = atomSet.getAtom(1);
        IAtom atom2 = atomSet.getAtom(2);
        dr12.Ev1Mv2(atom1.getPosition(),atom0.getPosition());
        dr23.Ev1Mv2(atom2.getPosition(),atom1.getPosition());
        boundary.nearestImage(dr12);
        boundary.nearestImage(dr23);
        double costheta = -dr12.dot(dr23)/Math.sqrt(dr12.squared()*dr23.squared());
        double dtheta;
        // machine precision can give us numbers with magnitudes slightly greater than 1
        if (costheta > 1) {
            dtheta = 0;
        }
        else if (costheta < -1) {
            dtheta = Math.PI;
        }
        else {
            dtheta = Math.acos(costheta);
        }
        dtheta -= angle;
        return 0.5*epsilon*dtheta*dtheta;
    }

    /**
     * Sets the nominal bond angle (in radians)
     */
    public void setAngle(double newAngle) {
        angle = newAngle;
    }
    
    /**
     * Returns the nominal bond angle (in radians)
     */
    public double getAngle() {
        return angle;
    }
    
    public Dimension getAngleDimension() {
        return Angle.DIMENSION;
    }

    /**
     * Sets the characteristic energy of the potential
     */
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    /**
     * Returns the characteristic energy of the potential
     */
    public double getEpsilon() {
        return epsilon;
    }
    
    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public Vector[] gradient(IAtomList atoms) {
        IAtom atom0 = atoms.getAtom(0);
        IAtom atom1 = atoms.getAtom(1);
        IAtom atom2 = atoms.getAtom(2);
        dr12.Ev1Mv2(atom1.getPosition(),atom0.getPosition());
        dr23.Ev1Mv2(atom2.getPosition(),atom1.getPosition());
        boundary.nearestImage(dr12);
        boundary.nearestImage(dr23);
        double dr12_23 = 1.0/Math.sqrt(dr12.squared()*dr23.squared());
        double costheta = -dr12.dot(dr23)*dr12_23;
        // machine precision can give us numbers with magnitudes slightly greater than 1
        if (costheta > 0.999 || costheta < -0.999) {
            // equation is 0/0
            gradient[0].E(0);
            gradient[1].E(0);
            gradient[2].E(0);
            return gradient; 
        }
        double dtheta = Math.acos(costheta) - angle;
        gradient[0].Ea1Tv1(dr12_23, dr23);
        gradient[0].PEa1Tv1(costheta/dr12.squared(), dr12);
        gradient[0].TE(-epsilon*dtheta/Math.sqrt(1.0-costheta*costheta));

        gradient[2].Ea1Tv1(-dr12_23, dr12);
        gradient[2].PEa1Tv1(-costheta/dr23.squared(), dr23);
        gradient[2].TE(-epsilon*dtheta/Math.sqrt(1.0-costheta*costheta));
        
        gradient[1].Ea1Tv1(-1, gradient[0]);
        gradient[1].PEa1Tv1(-1, gradient[2]);
        
        return gradient;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    protected final Vector dr12, dr23;
    protected Boundary boundary;
    protected double angle;
    protected double epsilon;
    private static final long serialVersionUID = 1L;
    protected final Vector[] gradient;
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        RandomNumberGenerator random = new RandomNumberGenerator();
        P3BondAngle potential = new P3BondAngle(space);
        potential.setEpsilon(1.0);
        potential.setAngle(Math.PI/2.0);
        Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
        RandomPositionSourceRectangular positionSource = new RandomPositionSourceRectangular(space, random);
        positionSource.setBox(box);
        potential.setBox(box);
        Atom atom0 = new Atom(space);
        atom0.getPosition().setX(0, 1);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        int n = 100;
        double oldU = 0;
        double oldoldU = 0;
        double U = 0;
        Vector oldGradient = space.makeVector();
        Vector gradient = space.makeVector();
        Vector dr = space.makeVector();
        for (int i=0; i<n+1; i++) {
            oldoldU = oldU;
            oldU = U;
            oldGradient.E(gradient);

            double theta = i*Math.PI/n;
            atom2.getPosition().setX(0, Math.cos(theta));
            atom2.getPosition().setX(1, Math.sin(theta));
            U = potential.energy(atoms);
            gradient.E(potential.gradient(atoms)[2]);
//            System.out.println(theta+" "+potential.energy(atoms)+" "+potential.gradient(atoms)[2]);
            dr.setX(0, -Math.cos((i-2)*Math.PI/n));
            dr.setX(1, -Math.sin((i-2)*Math.PI/n));
            dr.PE(atom2.getPosition());
            
            System.out.println((i-1)*Math.PI/n+" "+oldGradient.dot(dr)+" "+(U-oldoldU));
        }

        System.out.println("all random");
        for (int i=0; i<n; i++) {
            atom0.getPosition().E(positionSource.randomPosition());
            atom1.getPosition().E(positionSource.randomPosition());
            atom2.getPosition().E(positionSource.randomPosition());
            
            U = potential.energy(atoms);

            int iRand = random.nextInt(3);
            iRand = 1;
            IAtom atom = atoms.getAtom(iRand);
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
