package etomica.models.nitrogen;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
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
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Energy;
import etomica.units.Kelvin;
import etomica.util.FunctionGeneral;



/**
 * Lattice sum class for Alpha-phase Nitrogen
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicAlphaNitrogenModelLatticeSum extends Simulation{

	
	public HarmonicAlphaNitrogenModelLatticeSum(ISpace space, int numMolecule, double density, int iLayer) {
		super(space);
				
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
		double unitCellLength = Math.pow(numMolecule/density, 1.0/3.0)/nCell;//5.661;
//		System.out.println("a: " + unitCellLength);
//		System.out.println("nCell: " + nCell);

		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		SpeciesN2 species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		SpeciesN2B ghostSpecies = new SpeciesN2B(space);
		ghostSpecies.setConformation(conformation);
		addSpecies(ghostSpecies);
		
		
		IBox box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		IBox ghostBox = new Box(space);
		addBox(ghostBox);
		ghostBox.setNMolecules(ghostSpecies, 1);
		
		int [] nCells = new int[]{1,1,1};
		Primitive primitive = new PrimitiveCubic(space, nCell*unitCellLength);
	
		CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.initializeCoordinates(nCells);
	
		double rCScale = 0.475;
		double rC = 1000;//box.getBoundary().getBoxSize().getX(0)*rCScale;
		//System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		
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
		latticeSum.setMaxLatticeShell(iLayer);
		
		double sum = 0;
	    double basisDim = lattice.getBasis().getScaledCoordinates().length;
		DataGroupLSC data = (DataGroupLSC)latticeSum.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        double latEnergy = 0.5*sum/basisDim;
        System.out.println("lattice energy [sim unit]:  " + latEnergy + " ;[K]: " + Kelvin.UNIT.fromSim(latEnergy));
	
	}
	
	public static void main (String[] args){
		
		int numMolecule =4;
		double density = 0.025;

		System.out.println("Lattice Energy Calculation of Alpha-phase Nitrogen");
		System.out.println("Using lattice sum method with truncation of 1000A");
		System.out.println("with density of:" + density);
		
		for (int i=90; i<121; i++){	
			System.out.print(i+" ");
			HarmonicAlphaNitrogenModelLatticeSum test = new HarmonicAlphaNitrogenModelLatticeSum(Space3D.getInstance(3), numMolecule, density, i);
		}

	}
	
	private static final long serialVersionUID = 1L;
}
