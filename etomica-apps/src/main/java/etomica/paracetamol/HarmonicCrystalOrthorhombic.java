/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.WriteConfigurationDLPOLY;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.paracetamol.LatticeSumCrystalParacetamol.DataGroupLSCParacetamol;
import etomica.potential.P2DLPOLY;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Energy;
import etomica.units.Kelvin;
import etomica.util.Arrays;

/**
 * Properties of a system of monatomic molecules occupying a lattice and interacting according
 * to a spherically-symmetric pair potential.  Properties are given by a lattice-dynamics treatment.
 * 
 * @author Tai Tan
 *
 */
public class HarmonicCrystalOrthorhombic {

    public HarmonicCrystalOrthorhombic(int[] nCells, Primitive primitive,
    		     Basis basis, Box box,
    		     CoordinateDefinitionParacetamol coordinateDefinitionParacetamol,
    		     Space _space) {
        this.nCells = nCells.clone();
        this.coordinateDefinitionParacetamol = coordinateDefinitionParacetamol;
        this.space = _space;
        potential = new P2DLPOLY(space);
        lattice = new BravaisLatticeCrystal(primitive, basis);
        normalModes = new NormalModesPotentialParacetamol(nCells, primitive, basis, space);
        normalModes.setBox(box);
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    public double getLatticeEnergy() {
    	
    	/*
    	 * Lattice Energy is inner class that compute the pair energy
    	 *  between the two particular paracetamol molecules
    	 */
    	LatticeEnergy latticeEnergy = new LatticeEnergy(potential, coordinateDefinitionParacetamol);
    	
    	LatticeSumCrystalParacetamol summer = new LatticeSumCrystalParacetamol(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        summer.setK(lattice.getSpace().makeVector());
        summer.setCoordinateDefinitionParacetamol(coordinateDefinitionParacetamol);
//        aSystem.out.println("\n k:"+kVector.toString());
        
        double sum = 0;
        double basisDim = lattice.getBasis().getScaledCoordinates().length;
        
        DataGroupLSCParacetamol data = (DataGroupLSCParacetamol)summer.calculateSum(latticeEnergy);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        return 0.5*sum/basisDim;
//            a((Tensor)sum[0]).map(chopper);
//            a((Tensor)sum[1]).map(chopper);
//            a((Tensor)sum[0]).ME(sum0);
//            aSystem.out.println(sum[0].toString());
 //           aSystem.out.println();
//            aSystem.out.println(sum[1].toString());
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
        sumA += getLatticeEnergy(); //Get the lattic energy from the minimized configuration
        return sumA;
    }
    
    public void setCellDensity(double newDensity) {
        double oldVolume = lattice.getPrimitive().unitCell().getVolume();
        double scale = newDensity * oldVolume;
        Primitive primitive = lattice.getPrimitive();
        primitive.scaleSize(1.0/Math.pow(scale, 1.0/lattice.getSpace().D()));
//       a normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
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
        double temperature = Kelvin.UNIT.toSim(10);
        
//        double rho = 1.0;
//        int maxLatticeShell = 49;
//        aPrimitive primitive = new PrimitiveFcc(Space3D.getInstance());
//        aBasis basis = new BasisMonatomic(Space3D.getInstance());

        if (args.length > 0) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[0]));
        }
        
        Space sp = Space3D.getInstance();
        MCParacetamolOrthorhombicDLMULTI sim = new MCParacetamolOrthorhombicDLMULTI(sp,16,temperature, 3, new int[] {1,1,2});
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();
        
        
        //a double[] u = new double[]{-8.519546134513288, -6.357161053952702, -7.449621419035304, -0.003932863936184616, -0.0016776660890108416, -0.013522565626184084, -8.726086341054327, -5.722820924623122, -7.451563204983988, -0.005551030882875035, -0.008781239888660329, 0.005654659111204717, -8.737722525835544, -6.358584368530293, -7.312474753722753, -0.004632407803175042, 0.0073903160411125865, -0.0022154770629279168, -8.513251161366203, -5.718672928462283, -7.308964948201192, -0.008253140936357385, -0.0016033969986892181, -0.007589294260920792, -8.726675083680286, -5.711819873904577, -7.345829560494344, 0.0012238914448621102, 0.009773882446509393, 0.0029289771411575016, -8.4909577926205, -6.372620681124473, -7.328828834258954, 0.008259322829339406, 7.612944355468755E-4, 0.015416013389648568, -8.536258682509295, -5.721912267935009, -7.4175804780220504, 0.0048799658666174375, 0.01669560616620193, -0.004313176827043905, -8.741502278415917, -6.38040790144197, -7.441136801296486, 0.0043516618097694075, 0.010742961480577531, 0.0048663615734551415}; 
      
        BasisCell[] cell = sim.coordDef.getBasisCells() ;
      
        
        // for (int i=0; i<cell.length; i++){
        //  	coordinateDefinition.setToU(cell[i].molecules, u);
        // }
        
        WriteConfigurationDLPOLY configDLPOLY = new WriteConfigurationDLPOLY();
        configDLPOLY.setConfName("CONFIG");
        // configDLPOLY.setBox();
        configDLPOLY.getElementHash().put(HydrogenP.INSTANCE, "HP");
   
        //a final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);

        HarmonicCrystalOrthorhombic harmonicCrystal = new HarmonicCrystalOrthorhombic(new int[]{1,1,2}, sim.primitive, basis, sim.box, sim.coordDef, sp);
        //harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
        harmonicCrystal.normalModes.setCoordinateDefinition(sim.coordDef);
        harmonicCrystal.normalModes.setPotentialMaster(sim.potentialMaster);
        harmonicCrystal.setMaxLatticeShell(harmonicCrystal.maxLatticeShell);
   
        
        //System.out.println("Density: " + rho);
        System.out.println("Temperature: " + temperature);
        
        double uEnergy = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy(temperature);
        System.out.println("Lattice Energy: " + uEnergy);
        //double uEos = LennardJones.uStaticFcc(rho);
        //System.out.println("Lattice energy from EOS: " + uEos);
        System.out.println("Helmholtz: " + a);
        //double aEos = LennardJones.aResidualFcc(T,rho) + T*Math.log(rho) - 1.0*T;
        //System.out.println("Helmholtz from EOS: " + aEos);
        
//        a double latticeConstant = 1.0;
//        a primitive = new PrimitiveHexagonal(Space3D.getInstance(), latticeConstant, Math.sqrt(8.0/3.0)*latticeConstant);
//        a basis = new BasisHcp();
//        a harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential);
//        a harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
//        a harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
//        a u = harmonicCrystal.getLatticeEnergy();
//        a System.out.println("Lattice energy (HCP): "+u);
    }
    
    
    private NormalModesPotentialParacetamol normalModes;
    private BravaisLatticeCrystal lattice;
    private int[] nCells;
    private int maxLatticeShell;
    private final CoordinateDefinitionParacetamol coordinateDefinitionParacetamol;
    private final P2DLPOLY potential;
    private final Space space;
    private static final long serialVersionUID = 1L;
    
    public static final class LatticeEnergy implements LatticeEnergyParacetamol {
    	
    	public LatticeEnergy(P2DLPOLY potential, CoordinateDefinitionParacetamol coordinateDefinitionParacetamol){
    		this.potential = potential;
    		potential.getConfigP2DLPOLY().getElementHash().put(HydrogenP.INSTANCE, "HP");
    		potential.setBox(coordinateDefinitionParacetamol.getBox());
    		
    		this.coordinateDefinitionParacetamol = coordinateDefinitionParacetamol;
    		singleMolDim = coordinateDefinitionParacetamol.getCoordinateDim()
  		  	/coordinateDefinitionParacetamol.getBasisCells()[0].molecules.getMoleculeCount();
    		dataInfo = new DataInfoDouble("Lattice Energy", Energy.DIMENSION);
    		pairEnergy = new DataDouble();
    		tag = new DataTag();
   
    		
    	}
    
		public IData getData() {
			
			double[] u = new double[coordinateDefinitionParacetamol.getCoordinateDim()];
			
			IMoleculeList molecules = coordinateDefinitionParacetamol.getBasisCells()[0].molecules; 
			coordinateDefinitionParacetamol.setToU(molecules, u);
			
			MoleculePair pair = new MoleculePair(molecules.getMolecule(indexj), molecules.getMolecule(indexjp));
			
			pairEnergy.x = potential.energy(pair);
			
			return pairEnergy;
         }
		
		public void setMolecule(int indexj, int indexjp){
			this.indexj = indexj;
			this.indexjp = indexjp;
			
		}

		public IDataInfo getDataInfo() {
			 return dataInfo;
		}

		public DataTag getTag(){
			return tag;
		}
		
		private final CoordinateDefinitionParacetamol coordinateDefinitionParacetamol;
		private final P2DLPOLY potential;
		private final DataTag tag;
		private final DataDouble pairEnergy;
		private final DataInfo dataInfo;
		private final int singleMolDim;
		private int indexj, indexjp;
	}
    
    
    
}
