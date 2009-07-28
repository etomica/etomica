 package etomica.normalmode;

import etomica.data.DataInfo;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Energy;
import etomica.util.Arrays;
import etomica.util.FunctionGeneral;

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

    public HarmonicCrystalSoftSphereFCC(int[] nCells, Primitive primitive,
    		    Basis basis, Potential2SoftSpherical potential, ISpace _space, String file_Name) {
        this.potential = potential;
        this.nCells = (int[])nCells.clone();
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);
        normalModes = new NormalModesPotential(nCells, primitive, basis, potential, space);
        setFileName(file_Name);
        setMaxLatticeShell(49);
    }
    
    public double getLatticeEnergy() {
        FunctionGeneral function = new FunctionGeneral() {
            public IData f(Object obj) {
                data.x = potential.u(((Vector3D)obj).squared());
                return data;
            }
            public IDataInfo getDataInfo() {
                return dataInfo;
            }
            final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
            final DataDouble data = new DataDouble();
        };
        LatticeSumCrystal summer = new LatticeSumCrystal(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        summer.setK(lattice.getSpace().makeVector());
//        System.out.println("\n k:"+kVector.toString());
        double sum = 0;
        double basisDim = lattice.getBasis().getScaledCoordinates().length;
        DataGroupLSC data = (DataGroupLSC)summer.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        return 0.5*sum/basisDim;
//            ((Tensor)sum[0]).map(chopper);
//            ((Tensor)sum[1]).map(chopper);
//            ((Tensor)sum[0]).ME(sum0);
//            System.out.println(sum[0].toString());
 //           System.out.println();
//            System.out.println(sum[1].toString());
    }
    
    public double getHelmholtzFreeEnergy(double temperature) {
        
        int D = lattice.getSpace().D();
        int cellCount = 1;
        int differ = 1;
        for(int i=0; i<D; i++) {
            cellCount *= nCells[i];
            if(nCells[i] % 2 == 0) {
                differ *= 2;
            }
        }
        int basisDim = lattice.getBasis().getScaledCoordinates().length;
        int moleculeCount = cellCount*basisDim;
        System.out.println("Space dimension: " + D);
        System.out.println("cellCount: "+cellCount);
        System.out.println("basisDim: "+basisDim);
        System.out.println("moleculeCount: "+moleculeCount);
        double jacobian = 0.5*D*(basisDim*(cellCount - differ)*Math.log(2.0) - Math.log(cellCount));
        System.out.println("differ, jacobian: " + differ + "\t" + jacobian);
        
        normalModes.setFileName(getFileName());
        double[][] omega2 = normalModes.getOmegaSquared(null);//need to change signature of this method
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        
        System.out.println("coeffs: "+Arrays.toString(coeffs));
        double sumA = 0.0;
        double normalModeSum = 0.0;
        double omega2zeroCount = 0;
        for(int k=0; k<omega2.length; k++) {
            double coeff = coeffs[k];
            for(int i=0; i<omega2[k].length; i++) {
                if(omega2[k][i] > 1.e-9) {
                    sumA += coeff*Math.log(omega2[k][i]*coeff/(temperature*Math.PI));
                    normalModeSum += coeff;
                } else {
                    omega2zeroCount++;
                }
            }
        }

        System.out.println("omega2==0 count: "+omega2zeroCount);
        System.out.println("2*normalModeSum + D: " + (2*normalModeSum+D));
        System.out.println("D * moleculeCount: " + (D*moleculeCount));
        sumA -= jacobian;
        sumA /= moleculeCount;
        sumA *= temperature;
        sumA += getLatticeEnergy();
        return sumA;
    }
    
    public void setCellDensity(double newDensity) {
        double oldVolume = lattice.getPrimitive().unitCell().getVolume();
        double scale = newDensity * oldVolume;
        Primitive primitive = lattice.getPrimitive();
        primitive.scaleSize(1.0/Math.pow(scale, 1.0/lattice.getSpace().D()));
//        normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
        normalModes = new NormalModesPotential(nCells, primitive, lattice.getBasis(), potential, space);
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }
    
    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    public void setMaxLatticeShell(int maxLatticeShell) {
        this.maxLatticeShell = maxLatticeShell;
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }

	public String getFileName() {
		return fileName;
	}

	public void setFileName(String filename) {
		this.fileName = filename;
	}
    
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
        
        
        
        HarmonicCrystalSoftSphereFCC harmonicCrystal = new HarmonicCrystalSoftSphereFCC(nCells, primitive, basis, potential, sp, fileName);
       
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
        harmonicCrystal.setMaxLatticeShell(harmonicCrystal.maxLatticeShell);

        double u = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        
        System.out.println("\nLattice Energy: " + u);
        System.out.println("Helmholtz Free Energy at T"+temperature+ " is: "+a);
        System.out.println("Harmonic-reference free energy: "+ (a-u));
      
        System.out.println("\nCalcHarmonicA from file (Temperature-independent)");
        CalcHarmonicA calcHarmonicA = new CalcHarmonicA();
        calcHarmonicA.doit(fileName, 3, 1.0, temperature, basis.getScaledCoordinates().length, nC*nC*nC);
        
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
    
    
    private NormalModesPotential normalModes;
    private BravaisLatticeCrystal lattice;
    private int[] nCells;
    private int maxLatticeShell;
    private Potential2SoftSpherical potential;
    private final ISpace space;
    private String fileName;
    private static final long serialVersionUID = 1L;
    
}
