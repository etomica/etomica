package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.Degree;
import etomica.util.FunctionMultiDimensionalDifferentiable;
import etomica.util.numerical.ArrayReader1D;
import etomica.util.numerical.ConjugateGradientMultiDimensional;

/**
 * Lattice Energy Minimization Routine for beta-phase nitrogen using
 *  Conjugate Gradient. 
 *  
 * @author taitan
 *
 */
public class MinimizeBetaNitrogenLatticeParameterCG extends Simulation implements FunctionMultiDimensionalDifferentiable {
	
	public MinimizeBetaNitrogenLatticeParameterCG(ISpace space, int[] nC, double density){
		super(space);
		
		double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		int numMolecule = nC[0]*nC[1]*nC[2]*2;
		
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, nC);
		
		SpeciesN2 species = new SpeciesN2(space);
		addSpecies(species);
		
		IBox box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);
		
		IVector[] boxDim = new IVector[3];
		boxDim[0] = space.makeVector(new double[]{nC[0]*aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1]*aDim*Math.cos(Degree.UNIT.toSim(60)), nC[1]*aDim*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2]*cDim});
		
		Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
		Primitive primitive = new PrimitiveHexagonal(space, nC[0]*aDim, nC[2]*cDim);
		
		box.setBoundary(boundary);
		double rC = nC[0]*aDim*0.475;
		P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		potentialMaster = new PotentialMaster();
	    potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		coordinateDefinition = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDefinition.setIsBeta();
		coordinateDefinition.setOrientationVectorBeta(space);
		coordinateDefinition.initializeCoordinates(new int[]{1,1,1});
		
		derivativeFunc = new FiniteDifferenceDerivativeCGNitrogenBeta(box, potentialMaster, coordinateDefinition);
	}
	
	public double f(double[] u){
		return derivativeFunc.f(u);
	}
	
	public double df(int[] d, double[] u){
		return derivativeFunc.dphi_du(d, u);
	}
	
	public int getDimension(){
		return this.getDimension();
	}
	
	public static void main(String args[]){
		String filename = "/usr/users/taitan/workspace/etomica.Apps1/bin/inputd0.0240";
		double density = 0.0240; 
		int nCells = 8;
		
	    if(args.length > 0){
			filename = args[0];
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			nCells = Integer.parseInt(args[2]);
		}
		int[] nC = new int[]{nCells,nCells,nCells};
		int numMolecules = nC[0]*nC[1]*nC[2]*2;
		
		System.out.println("Running Lattice Energy Minimization for Beta-Phase Nitrogen");
		System.out.println("with density: "+density + " with numMolecules: "+numMolecules+"\n");
		
		double[] parameters = new double[20];
		double[][] paramFromFile = ArrayReader1D.getFromFile(filename);
		int k=0;
		for (int i=0; i<paramFromFile.length;i++){
			for (int j=0; j<paramFromFile[0].length;j++){
				
				parameters[k]=paramFromFile[i][j];
				k++;
			}	
		}
		
//		double[] parameters = new double[]{
//				-1.4621203267796472E-4, -7.171444595461441E-5, -2.37287496966035E-5, -2.910494246533862E-5, -2.1177506743927406E-5, 
//				-1.4659318210114304E-4, -9.306680203123414E-5, -2.9620966817594917E-5, -5.5768512729374026E-5, -5.2196511228923316E-5, 
//				-1.3068268031769714E-4, -6.723384250248303E-5, -2.3977416456710433E-5, 2.0756040199971303E-4, -1.7564319690730054E-4, 
//				-1.377420510336E-4, -1.12939367499009E-4, -2.373467463205741E-5, 1.0206286755957938E-5, -1.8563495206484917E-5
//		};
//		double[] parameters = new double[]{
//				0.13472369149746394, -7.777011074268426E-5, 2.2896362375169623E-4, -3.519301914033887E-5, -2.056152153036304E-5, 
//				0.13471923410456724, -1.1020893715247124E-4, 2.227170296935758E-4, -1.0088482846707517E-4, -7.004772945742876E-5, 
//				0.13474085174082737, -7.362751408992318E-5, 2.2355286629845394E-4, 2.4940426703535894E-4, -1.870678763479276E-4, 
//				0.1347235701316199, -1.1600382173691587E-4, 2.3020012530647829E-4, 1.643752779412067E-5, 3.0229264499046215E-6
//		};
				
		MinimizeBetaNitrogenLatticeParameterCG testFunction = new MinimizeBetaNitrogenLatticeParameterCG(Space.getInstance(3), nC, density);
		System.out.println("Initial Energy value: "+ testFunction.f(parameters));
	
		ConjugateGradientMultiDimensional conjugateGradient = new ConjugateGradientMultiDimensional();
		conjugateGradient.conjugateGradient(parameters, 1.2e-9, testFunction);
		
		double minLatEnergy =conjugateGradient.getFunctionMinimimumValue();
		System.out.println("Minimum Function value is: "+ minLatEnergy+" "+minLatEnergy/numMolecules);
		
		parameters = conjugateGradient.getMinimumCoordinates();
//		for(int i=0; i<parameters.length; i++){
//			System.out.print(parameters[i]+", ");
//			if (i>0 && (i+1)%5==0)System.out.println();
//		}

		try {
			FileWriter fileWriter = new FileWriter(filename,false);
			
			for (int i=0; i<parameters.length; i++){
				fileWriter.write(parameters[i]+" ");
				
				if(i>0&&(i+1)%5==0){
					fileWriter.write("\n");
						
				}
			}
			//fileWriter.write(minLatEnergy/numMolecules+"\n");
			fileWriter.close();
			
		} catch(IOException e){
			throw new RuntimeException("Failed to write file" + e);
		
		}
	}
	
	protected FiniteDifferenceDerivativeCGNitrogenBeta derivativeFunc;
	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected MeterPotentialEnergy meterPotential;
	protected PotentialMaster potentialMaster;
	
	private static final long serialVersionUID = 1L;

}
