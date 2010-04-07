 package etomica.normalmode;

import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
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
 * HCP Crystal Structure
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicCrystalSoftSphereHCP {

    public static void main(String[] args) {
        double rho = 1.1964;
        int exponent = 12;
        int maxLatticeShell = 6;
        int nC =4;
        int numAtom = nC*nC*nC*2;
        
        double temperature = 1.0;
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
        
        double a = Math.pow(Math.sqrt(2)/rho, 1.0/3.0);
        double c = Math.sqrt(8.0/3.0)*a;
                
        Primitive primitive = new PrimitiveHexagonal(Space3D.getInstance(),nC*a,nC*c);
        Basis basisHCP = new BasisHcp();
        
        Basis basis = new BasisBigCell(Space.getInstance(3), basisHCP, new int[]{nC,nC,nC});
        
        ISpace sp = Space3D.getInstance();
        final Potential2SoftSpherical potential = new P2SoftSphere(sp, 1.0, 1.0, exponent);
        double rc = (maxLatticeShell*2)*nC*Math.pow(4.0/rho, 1.0/3.0)*0.495;
        System.out.println("truncation at "+rc);
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(sp, potential, rc);
        
        int[] nCells = new int[] {1, 1, 1};
        long startTime = System.currentTimeMillis();
        
        // we want an untruncated lattice energy, so extend lattice sum by 2 more cells
        double rcU = ((maxLatticeShell+2)*2)*nC*Math.pow(4.0/rho, 1.0/3.0)*0.495;
        pTruncated.setTruncationRadius(rcU);
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential, sp);
        
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);

        harmonicCrystal.setMaxLatticeShell(maxLatticeShell+2);
        double u = harmonicCrystal.getLatticeEnergy();
        System.out.println("lattice energy: " + u);
        System.exit(1);

        // now actually calculate the free energy
        harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, pTruncated, sp);
       
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);

        harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
        harmonicCrystal.getNormalModes().setFileName(fileName);
        
        System.out.println("0.001  " + harmonicCrystal.getHelmholtzFreeEnergy(0.001));
        System.out.println("1.2147 " + harmonicCrystal.getHelmholtzFreeEnergy(1.2147));
        for (int i=1; i<17; i++){
        	System.out.println((i*0.1)+" " + harmonicCrystal.getHelmholtzFreeEnergy(i*0.1));
        }
        
        double f = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        
        System.out.println("\nLattice Energy: " + u);
        System.out.println("Helmholtz Free Energy at T"+temperature+ " is: "+f);
        System.out.println("Harmonic-reference free energy: "+ (f-u));
      
        System.out.println("\nCalcHarmonicA from file (Temperature-independent)");
        CalcHarmonicA.doit(harmonicCrystal.getNormalModes(), 3, temperature, basis.getScaledCoordinates().length);

//        long endTime = System.currentTimeMillis();
//        System.out.println("End Time: " + endTime);
//        System.out.println("Time taken: " + (endTime - startTime));
        

    }
}
