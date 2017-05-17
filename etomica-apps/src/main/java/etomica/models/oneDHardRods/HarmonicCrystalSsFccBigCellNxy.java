/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.CalcHarmonicA;
import etomica.normalmode.HarmonicCrystal;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Properties of a system of monatomic molecules occupying a lattice and interacting according
 * to a soft-sphere spherically-symmetric pair potential.  Properties are given by a lattice-dynamics treatment.
 * 
 * FCC Crystal Structure
 * 
 * @author kofke & Tai Tan
 *
 */
public class HarmonicCrystalSsFccBigCellNxy {

    public static void main(String[] args) {
        double rho = 1.1964;
        int exponent = 12;
        int maxLatticeShell = 1;
        double primitiveLength = 1.4952983783596188;
        System.out.println("MaxLatticeShell is "+maxLatticeShell);
        int nC = 8;
        int[] shape = new int[] {nC, nC, nC};
        int numAtom = 1;
        for(int i = 0; i < 3; i++){
            numAtom *= shape[i];
        }
        numAtom *= 4;
        
        double temperature = 0.01;
        String fileName = "inputSSDB_BC"+ numAtom;
        
        if (args.length > 0) {
            rho = Double.parseDouble(args[0]);
        }
        if (args.length > 1) {
            exponent = Integer.parseInt(args[1]);
        }
        if (args.length > 2) {
            nC = Integer.parseInt(args[2]);
        }
        if (args.length > 3) {
            temperature = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
            fileName = args[4];
        }
        
        System.out.println("FCC Harmonic Crystal with n = " + exponent);
        System.out.println("numAtom: " + numAtom);
        
        double[] lengths = new double[3];
        lengths[0] = shape[0] * primitiveLength;
        lengths[1] = shape[1] * primitiveLength;
        lengths[2] = shape[2] * primitiveLength;
        
        Primitive primitive = new PrimitiveOrthorhombic(Space3D.getInstance(), 
                lengths[0], lengths[1], lengths[2] );
        Basis basisFCC = new BasisCubicFcc();
        
        Basis basis = new BasisBigCell(Space.getInstance(3), basisFCC, shape);
        
        Space sp = Space3D.getInstance();
        final Potential2SoftSpherical potential = new P2SoftSphere(sp, 1.0, 1.0, exponent);
        double rc = 15 ; //(maxLatticeShell*2)*nC*Math.pow(4.0/rho, 1.0/3.0)*0.495;
        //nan  Set for myself
        rc = 0.495 * 2 * primitiveLength;
        if( numAtom > 255) { rc = 2.2;}
        System.out.println("truncation "+rc);
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(sp, potential, rc);
        
        int[] nCells = new int[] {1, 1, 1};
      
        pTruncated.setTruncationRadius(rc);
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(rho, nCells, primitive, basis, pTruncated, sp);
        harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
        harmonicCrystal.getNormalModes().setFileName(fileName);
        
        double u = harmonicCrystal.getLatticeEnergy();
        System.out.println("lattice energy: " + u);
    
        double a = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        
        System.out.println("\nLattice Energy: " + u);
        System.out.println("Helmholtz Free Energy at T"+temperature+ " is: "+a);
        System.out.println("Harmonic-reference free energy: "+ (a-u));
      
        System.out.println("\nCalcHarmonicA from file (Temperature-independent)");
        CalcHarmonicA.doit(harmonicCrystal.getNormalModes(), 3, temperature, basis.getScaledCoordinates().length);

    }
}
