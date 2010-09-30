package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.DataInfo;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.models.nitrogen.LatticeSumCrystalMolecular.DataGroupLSC;
import etomica.normalmode.BasisBigCell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Energy;
import etomica.util.FunctionGeneral;



/**
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicAlphaNitrogenModelLatticeSum extends Simulation{

	
	public HarmonicAlphaNitrogenModelLatticeSum(ISpace space, int numMolecule, double density) {
		super(space);
		this.space = space;
		
		
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
		double unitCellLength = Math.pow(numMolecule/density, 1.0/3.0)/nCell;//5.661;
		System.out.println("a: " + unitCellLength);
		System.out.println("nCell: " + nCell);
		
		potentialMaster = new PotentialMaster();
				
		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		SpeciesN2 species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		SpeciesN2B ghostSpecies = new SpeciesN2B(space);
		ghostSpecies.setConformation(conformation);
		addSpecies(ghostSpecies);
		
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		IBox ghostBox = new Box(space);
		addBox(ghostBox);
		ghostBox.setNMolecules(ghostSpecies, 1);
		
		int [] nCells = new int[]{1,1,1};
	//	Boundary boundary = new BoundaryDeformablePeriodic(space,nCell*unitCellLength);
		Primitive primitive = new PrimitiveCubic(space, nCell*unitCellLength);
	
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.initializeCoordinates(nCells);
		
	
	//	box.setBoundary(boundary);
		double rCScale = 0.475;
		double rC =11;//box.getBoundary().getBoxSize().getX(0)*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		
		final P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		potential.setEnablePBC(false);
		
		FunctionGeneral function = new FunctionGeneral() {
			public IData f(Object obj) {
				data.x = potential.energy((IMoleculeList)obj);
				return data;
			}
			public IDataInfo getDataInfo() {
				return dataInfo;
			}
			final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
			final DataDouble data = new DataDouble();
		};
		
		BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basis);
		LatticeSumCrystalMolecular latticeSum = new LatticeSumCrystalMolecular(lattice, coordinateDef, ghostBox);
		latticeSum.setMaxLatticeShell(4);
		
		double sum = 0;
	    double basisDim = lattice.getBasis().getScaledCoordinates().length;
		DataGroupLSC data = (DataGroupLSC)latticeSum.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        System.out.println("lattice:  " + 0.5*sum/basisDim);
		System.exit(1);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

	}
	
	public static void main (String[] args){
		
		int numMolecule =4;
		double density = 0.025;
		HarmonicAlphaNitrogenModelLatticeSum test = new HarmonicAlphaNitrogenModelLatticeSum(Space3D.getInstance(3), numMolecule, density);

		CalcNumerical2ndDerivative cm2ndD = new CalcNumerical2ndDerivative(test.box, test.potentialMaster,test.coordinateDef);
		
		double[] newU = new double[test.coordinateDef.getCoordinateDim()];
		
		String fname = new String (numMolecule+"_2ndDer");
		try {
			FileWriter fileWriter = new FileWriter(fname);
			
			double value = 0;
			for (int i=0; i<newU.length; i++){
				for (int j=0; j<newU.length; j++){
					value = cm2ndD.d2phi_du2(new int[]{i,j}, newU);
					
					if(Math.abs(value) < 1e-6){
						value = 0.0;
					}
					fileWriter.write(value+ " ");
				}
				fileWriter.write("\n");
			}
			fileWriter.close();
			
		} catch (IOException e) {
			
		}
	
//		System.out.println("d2phi_du2: " + cm2ndD.d2phi_du2(new int[]{5,54}, newU));
//		System.out.println("d2phi_du2: " + cm2ndD.d2phi_du2(new int[]{54,5}, newU));

	}
	
	
	protected Box box;
	protected ISpace space;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected PotentialMaster potentialMaster;
	private static final long serialVersionUID = 1L;
}
