/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataInfo;
import etomica.data.FunctionData;
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
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.statmech.LennardJones;
import etomica.units.dimensions.Energy;

/**
 * Properties of a system of monatomic molecules occupying a lattice and interacting according
 * to a spherically-symmetric pair potential.  Properties are given by a lattice-dynamics treatment.
 * 
 * @author kofke
 *
 */
public class HarmonicCrystal {

    public HarmonicCrystal(double rho, int[] nCells, Primitive primitive, Basis basis, Potential2SoftSpherical potential, Space _space) {
        this.potential = potential;
        this.nCells = nCells.clone();
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);
        setCellDensity(rho/basis.getScaledCoordinates().length);
        normalModes = new NormalModesPotential(nCells, primitive, basis, potential, space);
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
//        System.out.println("Space dimension: " + D);
//        System.out.println("cellCount: "+cellCount);
//        System.out.println("basisDim: "+basisDim);
//        System.out.println("moleculeCount: "+moleculeCount);
        double Acom = 0.5*D*Math.log(moleculeCount);
//        System.out.println("differ, COM: " + differ + "\t" + Acom);

        double[][] omega2 = normalModes.getOmegaSquared();
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
//        System.out.println("coeffs: "+Arrays.toString(coeffs));
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

//        System.out.println("omega2==0 count: "+omega2zeroCount);
//        System.out.println("2*normalModeSum + D: " + (2*normalModeSum+D));
//        System.out.println("D * moleculeCount: " + (D*moleculeCount));
        sumA -= Acom;
        sumA /= moleculeCount;
        sumA *= temperature;
        //sumA += getLatticeEnergy();
        return sumA;
    }
    
    /**
     * This sets the primitive to be the size that will give us the desired
     * cell density.
     */
    private void setCellDensity(double newDensity) {
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
        double T = 1.0;
        double rho = 1.256;
//        Primitive primitive = new PrimitiveFcc(Space3D.getInstance());
//        Basis basis = new BasisMonatomic(Space3D.getInstance());

        Space sp = Space3D.getInstance();
        Primitive primitive = new PrimitiveCubic(sp);
        Basis basis = new BasisCubicFcc();
        
        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);

        int nC = 2;
        int[] nCells = new int[] {nC, nC, nC};
        
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(rho, nCells, primitive, basis, potential, sp);
        harmonicCrystal.setMaxLatticeShell(2);
        
        System.out.println("Density: " + rho);
        System.out.println("Temperature: " + T);
        
        double u = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy(T);
        System.out.println("Lattice Energy: " + u);
        double uEos = LennardJones.uStaticFcc(rho);
        System.out.println("Lattice energy from EOS: " + uEos);
        System.out.println("Helmholtz: " + a);
        double aEos = LennardJones.aResidualFcc(T,rho) + T*Math.log(rho) - 1.0*T;
        System.out.println("Helmholtz from EOS: " + aEos);
        
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
    private final Space space;
}
