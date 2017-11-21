/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.WriteConfigurationDLPOLY;
import etomica.box.Box;
import etomica.data.DataInfo;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.PotentialDLPOLY;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.units.dimensions.Energy;

import java.util.Arrays;

/**
 * Properties of a system of monatomic molecules occupying a lattice and interacting according
 * to a spherically-symmetric pair potential.  Properties are given by a lattice-dynamics treatment.
 * 
 * @author kofke
 *
 */
public class HarmonicCrystalMonoclinic {

    public HarmonicCrystalMonoclinic(int[] nCells, Primitive primitive, Basis basis, Box box, Space _space) {
        //this.potential = potential;
        this.nCells = nCells.clone();
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);
        normalModes = new NormalModesPotentialParacetamol(nCells, primitive, basis, space);
        normalModes.setBox(box);
        //setMaxLatticeShell(49);
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    public double getLatticeEnergy() {
        FunctionData<Object> function = new FunctionData<Object>() {
            public IData f(Object obj) {
                //data.x = potential.u(((Vector3D)obj).squared());
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
  ///////////////////////////////////////////////////////////////////////////////////////////////////////  
    
    
    
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
        
        double[][] omega2 = normalModes.getOmegaSquared();
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
        //sumA += getLatticeEnergy(); //Get the lattic energy from the minimized configuration
        return sumA;
    }
    
    public void setCellDensity(double newDensity) {
        double oldVolume = lattice.getPrimitive().unitCell().getVolume();
        double scale = newDensity * oldVolume;
        Primitive primitive = lattice.getPrimitive();
        primitive.scaleSize(1.0/Math.pow(scale, 1.0/lattice.getSpace().D()));
//        normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
        normalModes = new NormalModesPotentialParacetamol(nCells, primitive, lattice.getBasis(), space);
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }
    
    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    public void setMaxLatticeShell(int maxLatticeShell) {
        this.maxLatticeShell = maxLatticeShell;
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }

    public static void main(String[] args) {
        double temperature = Kelvin.UNIT.toSim(123);
        
//        double rho = 1.0;
//        int maxLatticeShell = 49;
//        Primitive primitive = new PrimitiveFcc(Space3D.getInstance());
//        Basis basis = new BasisMonatomic(Space3D.getInstance());

        if (args.length > 0) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[0]));
        }
        
        etomica.space.Space sp = Space3D.getInstance();
        MCParacetamolMonoclinicDLMULTI sim = new MCParacetamolMonoclinicDLMULTI(sp,32,temperature, 2, new int[] {2,2,2});
        BasisMonoclinicParacetamol basis = new BasisMonoclinicParacetamol();
      
        CoordinateDefinitionParacetamol coordinateDefinition = new CoordinateDefinitionParacetamol(sim, sim.box, sim.primitive, basis, sp);
        coordinateDefinition.setBasisMonoclinic();
        //coordinateDefinition.setConfiguration(configFile);
        
        int[] nCells = new int[] {2, 2, 2};
        coordinateDefinition.initializeCoordinates(nCells);
        double[] u = new double[]{0.19923028993600184, -0.10831241063851138, -0.3374206766423242, -0.056318915244514295, -0.08373011094015517, -0.20967425215989952, -0.15036406963107662, -0.06903390627642114, 0.2911023015207981, -0.062482609873595024, -0.0836259970912052, -0.17490727322325597, -0.19043886713958358, 0.09585326949021145, 0.29982577023219825, -0.052978731871725596, 0.0794450448846585, 0.20078353728718995, 0.1791946947288432, 0.07473598884040555, -0.33875779999181965, -0.06542404448301108, 0.0856334706980024, 0.20954733660556393};
        BasisCell[] cell = coordinateDefinition.getBasisCells() ;
        for (int i=0; i<cell.length; i++){
        	coordinateDefinition.setToU(cell[i].molecules, u);
        }
        
        WriteConfigurationDLPOLY configDLPOLY = new WriteConfigurationDLPOLY();
        configDLPOLY.setConfName("CONFIG");
        //configDLPOLY.setBox();
        configDLPOLY.getElementHash().put(HydrogenP.INSTANCE, "HP");
        
        PotentialDLPOLY potentialDLPOLY = new PotentialDLPOLY(Space3D.getInstance());
        potentialDLPOLY.setConfigDLPOLY(configDLPOLY);
        
        //final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);

        
        
        
        HarmonicCrystalMonoclinic harmonicCrystal = new HarmonicCrystalMonoclinic(nCells, sim.primitive, basis, sim.box, sp);
        //harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
        harmonicCrystal.normalModes.setCoordinateDefinition(coordinateDefinition);
        harmonicCrystal.normalModes.setPotentialMaster(sim.potentialMaster);
        harmonicCrystal.setMaxLatticeShell(harmonicCrystal.maxLatticeShell);
   
        
        //System.out.println("Density: " + rho);
        System.out.println("Temperature: " + temperature);
        
        //double uEnergy = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        //System.out.println("Lattice Energy: " + u);
        //double uEos = LennardJones.uStaticFcc(rho);
        //System.out.println("Lattice energy from EOS: " + uEos);
        System.out.println("Helmholtz: " + a);
        //double aEos = LennardJones.aResidualFcc(T,rho) + T*Math.log(rho) - 1.0*T;
        //System.out.println("Helmholtz from EOS: " + aEos);
        
//        double latticeConstant = 1.0;
//        primitive = new PrimitiveHexagonal(Space3D.getInstance(), latticeConstant, Math.sqrt(8.0/3.0)*latticeConstant);
//        basis = new BasisHcp();
//        harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential);
//        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
//        harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
//        u = harmonicCrystal.getLatticeEnergy();
//        System.out.println("Lattice energy (HCP): "+u);
    }
    
    
    private NormalModesPotentialParacetamol normalModes;
    private BravaisLatticeCrystal lattice;
    private int[] nCells;
    private int maxLatticeShell;
    private final Space space;
    private static final long serialVersionUID = 1L;
    
}
