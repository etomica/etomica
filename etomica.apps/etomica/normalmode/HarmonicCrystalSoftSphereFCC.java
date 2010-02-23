 package etomica.normalmode;

import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.ISpace;
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
public class HarmonicCrystalSoftSphereFCC {

    public static void main(String[] args) {
        double rho = 1.256;
        int exponent = 12;
        int maxLatticeShell = 2;
        int nC = 5;
        int numAtom = nC*nC*nC*4;
        
        double temperature = 0.001;
        String fileName = "inputSSDB"+ numAtom;
//        Primitive primitive = new PrimitiveFcc(Space3D.getInstance());
//        Basis basis = new BasisMonatomic(Space3D.getInstance());
        
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
        
        Primitive primitive = new PrimitiveCubic(Space3D.getInstance());
        Basis basisFCC = new BasisCubicFcc();
        
        Basis basis = new BasisBigCell(Space.getInstance(3), basisFCC, new int[]{nC,nC,nC});
        
        ISpace sp = Space3D.getInstance();
        final Potential2SoftSpherical potential = new P2SoftSphere(sp, 1.0, 1.0, exponent);
        double rc = (maxLatticeShell*2)*nC*Math.pow(4.0/rho, 1.0/3.0)*0.495;
        System.out.println("truncation at "+rc);
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(sp, potential, rc);
        
        int[] nCells = new int[] {1, 1, 1};
        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
        
        System.out.println("Running lattice-dynamics derivatives-based FCC soft-sphere simulation");
        System.out.println("Temperature: " + temperature);
        System.out.println("Density: " + rho);
        System.out.println("Exponent: " + exponent +"\n");
        
        // we want an untruncated lattice energy, so extend lattice sum by 2 more cells
        double rcU = ((maxLatticeShell+2)*2)*nC*Math.pow(4.0/rho, 1.0/3.0)*0.495;
        pTruncated.setTruncationRadius(rcU);
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, pTruncated, sp);
        
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);

        harmonicCrystal.setMaxLatticeShell(maxLatticeShell+2);
        double u = harmonicCrystal.getLatticeEnergy();


        // now actually calculate the free energy
        harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, pTruncated, sp);
       
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);

        harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
        harmonicCrystal.getNormalModes().setFileName(fileName);
        double a = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        
        System.out.println("\nLattice Energy: " + u);
        System.out.println("Helmholtz Free Energy at T"+temperature+ " is: "+a);
        System.out.println("Harmonic-reference free energy: "+ (a-u));
      
        System.out.println("\nCalcHarmonicA from file (Temperature-independent)");
        CalcHarmonicA.doit(harmonicCrystal.getNormalModes(), 3, temperature, basis.getScaledCoordinates().length);

        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken: " + (endTime - startTime));
        

    }
}
