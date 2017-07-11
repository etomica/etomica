/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Energy;
import etomica.data.FunctionData;
import etomica.util.ParameterBase;

/**
 * Properties of a system of monatomic molecules occupying a lattice and 
 * interacting according to a spherically-symmetric pair potential.  Properties 
 * are given by a lattice-dynamics treatment.
 * 
 * set up to use truncated soft sphere potential and output files 
 * 
 * @author cribbin off HarmonicCrystal of kofke
 * 
 */
public class HarmonicCrystalSsFccNxy {

    public HarmonicCrystalSsFccNxy(int[] nCells, Primitive primitive, Basis 
            basis, Potential2SoftSpherical potential, Space _space) {
        this.potential = potential;
        this.nCells = nCells.clone();
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);
        normalModes = new NormalModesPotential(nCells, primitive, basis, 
                potential, space);
    }

    public NormalModesPotential getNormalModes() {
        return normalModes;
    }

    public double getLatticeEnergy() {
        FunctionData<Object> function = new FunctionData<Object>() {
            public IData f(Object obj) {
                data.x = potential.u(((Vector3D)obj).squared());
                return data;
            }
            public IDataInfo getDataInfo() {
                return dataInfo;
            }
            final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice " +
            		"energy", Energy.DIMENSION);
            final DataDouble data = new DataDouble();
        };
        LatticeSumCrystal summer = new LatticeSumCrystal(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        summer.setK(lattice.getSpace().makeVector());
        double sum = 0;
        double basisDim = lattice.getBasis().getScaledCoordinates().length;
        DataGroupLSC data = (DataGroupLSC)summer.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        return 0.5*sum/basisDim;
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
        double Acom = 0.5*D*Math.log(moleculeCount);

        double[][] omega2 = normalModes.getOmegaSquared();
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        double sumA = 0.0;
        double normalModeSum = 0.0;
        double omega2zeroCount = 0;
        for(int k=0; k<omega2.length; k++) {
            double coeff = coeffs[k];
            for(int i=0; i<omega2[k].length; i++) {
                if(!Double.isInfinite(omega2[k][i])) {
                    sumA += coeff*Math.log(omega2[k][i]/(2*temperature*Math.PI));
                    normalModeSum += coeff;
                } else {
                    omega2zeroCount++;
                }
            }
        }

        sumA -= Acom;
        sumA /= moleculeCount;
        sumA *= temperature;
        return sumA;
    }
    
    public void setCellDensity(double newDensity) {
        double oldVolume = lattice.getPrimitive().unitCell().getVolume();
        double scale = newDensity * oldVolume;
        Primitive primitive = lattice.getPrimitive();
        primitive.scaleSize(1.0/Math.pow(scale, 1.0/lattice.getSpace().D()));
    }
    
    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    public void setMaxLatticeShell(int maxLatticeShell) {
        this.maxLatticeShell = maxLatticeShell;
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }

    public static void main(String[] args) {
        Params params = new Params();
        
        double T = params.T;
        double rho = params.rho;
        int[] nCells = params.shape;
        int nA = (nCells[0] * nCells[1] * nCells[2] * 4);
        String filename = params.filename + nA;
        double rc = params.rc ;
        
        Space sp = Space3D.getInstance();
        Potential2SoftSpherical potentialBase = new P2SoftSphere(sp, 1.0,
                1.0, 12);
        P2SoftSphericalTruncated potential = new P2SoftSphericalTruncated(
                sp, potentialBase, rc);
        potential.setTruncationRadius(rc);
        
        Primitive primitive = new PrimitiveCubic(sp);
        Basis basis = new BasisCubicFcc();
        
        double L = Math.pow(basis.getScaledCoordinates().length / rho, 1.0 / 3.0);
        primitive.scaleSize(L);
        
        HarmonicCrystalSsFccNxy harmonicCrystal = new HarmonicCrystalSsFccNxy(
                nCells, primitive, basis, potential, sp);
        harmonicCrystal.setMaxLatticeShell(3);
        harmonicCrystal.getNormalModes().setFileName(filename);
        
        System.out.println("Density: " + rho);
        System.out.println("Temperature: " + T);
        double a = harmonicCrystal.getHelmholtzFreeEnergy(T);
        double u = harmonicCrystal.getLatticeEnergy();
                
        System.out.println("\nLattice Energy: " + u);
        System.out.println("Helmholtz Free Energy at T "+T+ " is: "+a);
        System.out.println("Harmonic-reference free energy: "+ (a-u));
      
        System.out.println("\nCalcHarmonicA from file (Temperature-independent)");
        CalcHarmonicA.doit(harmonicCrystal.getNormalModes(), 3, T, nA);

    }
    
    private NormalModesPotential normalModes;
    private BravaisLatticeCrystal lattice;
    private int[] nCells;
    private int maxLatticeShell;
    private Potential2SoftSpherical potential;
    private final Space space;
    private static final long serialVersionUID = 1L;
    
    
    public static class Params extends ParameterBase {
        public double T = 0.01;
        public double rho = 1.1964;
        public int[] shape = new int[] {2, 2, 2};
        public String filename = "inputSSDB_";
        public double rc = 1.4803453945760225 ;
    }
    
}
