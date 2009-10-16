 package etomica.normalmode;

import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.ISpace;
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
        double rho = 1.116;
        int exponent = 12;
        int maxLatticeShell = 49;
        int nC =2;
        double temperature = 1.0;
        String fileName = "DB_FCC_n12_T16";
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
        Basis basis = new BasisCubicFcc();
        
        ISpace sp = Space3D.getInstance();
        final Potential2SoftSpherical potential = new P2SoftSphere(sp, 1.0, 1.0, exponent);

        int[] nCells = new int[] {nC, nC, nC};
        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
        
        System.out.println("Running lattice-dynamics derivatives-based FCC soft-sphere simulation");
        System.out.println("Temperature: " + temperature);
        System.out.println("Density: " + rho);
        System.out.println("Exponent: " + exponent +"\n");
        
        
        
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential, sp);
       
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);

        double u = harmonicCrystal.getLatticeEnergy();
        harmonicCrystal.getNormalModes().setFileName(fileName);
        double a = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        
        System.out.println("\nLattice Energy: " + u);
        System.out.println("Helmholtz Free Energy at T"+temperature+ " is: "+a);
        System.out.println("Harmonic-reference free energy: "+ (a-u));
      
        System.out.println("\nCalcHarmonicA from file (Temperature-independent)");
        CalcHarmonicA.doit(harmonicCrystal.getNormalModes(), 3, temperature, basis.getScaledCoordinates().length, nC*nC*nC);

        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken: " + (endTime - startTime));
        
//        double latticeConstant = 1.0;
//        primitive = new PrimitiveHexagonal(Space3D.getInstance(), latticeConstant, Math.sqrt(8.0/3.0)*latticeConstant);
//        basis = new BasisHcp();
//        harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential);
//        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
//        harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
//        u = harmonicCrystal.getLatticeEnergy();
//        System.out.println("Lattice energy (HCP): "+u);
    }
}
