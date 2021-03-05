/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.DataInfo;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.models.nitrogen.LatticeSumCrystalMolecular.DataGroupLSC;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.BasisBigCell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Joule;
import etomica.units.dimensions.Energy;


/**
 * Lattice sum class for Alpha-phase Nitrogen
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicAlphaNitrogenModelLatticeSum extends Simulation{

	
	public HarmonicAlphaNitrogenModelLatticeSum(Space space, int numMolecule, double density, double rC) {
		super(space);
				
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
		double unitCellLength = Math.pow(numMolecule/density, 1.0/3.0)/nCell;//5.661;
		System.out.println("a: " + unitCellLength);
//		System.out.println("nCell: " + nCell);

		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});
		
		SpeciesGeneral species = SpeciesN2.create(false);
		addSpecies(species);
		
		SpeciesGeneral ghostSpecies = SpeciesN2.create(false);
		addSpecies(ghostSpecies);
		
		
		Box box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);

		Box ghostBox = new Box(space);
		addBox(ghostBox);
		ghostBox.setNMolecules(ghostSpecies, 1);

		int[] nCells = new int[]{1, 1, 1};
		Primitive primitive = new PrimitiveCubic(space, nCell * unitCellLength);

		CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(getSpeciesManager(), box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.setIsDoLatticeSum();
		coordinateDef.initializeCoordinates(nCells);

		double rCScale = 0.475;
		//double rC = 1000;//box.getBoundary().getBoxSize().getX(0)*rCScale;
		//System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);

		final P2Nitrogen potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		potential.setEnablePBC(false);
		
		FunctionData<Object> function = new FunctionData<Object>() {
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
		
		double rX = unitCellLength*nCell;
		int nLayer = (int)Math.round(rC/rX + 0.5);
		System.out.println("rX: " + rX);
		System.out.println("nLayer: " + nLayer);
		
		BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basis);
		LatticeSumCrystalMolecular latticeSum = new LatticeSumCrystalMolecular(lattice, coordinateDef, ghostBox);
		latticeSum.setMaxLatticeShell(nLayer);
		
		double sum = 0;
	    double basisDim = lattice.getBasis().getScaledCoordinates().length;
		DataGroupLSC data = (DataGroupLSC)latticeSum.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        double latEnergy = 0.5*sum/basisDim;
        double avogradoConst = 6.0221415e23;
        System.out.println(density + " lattice energy [sim unit]:  " + latEnergy + " ;[kJ/mol]: " + Joule.UNIT.fromSim(latEnergy)*avogradoConst/1000);
	
	}
	
	public static void main (String[] args){
		
		int numMolecule =4;
		double density = 0.0230;
		double rC = 1200;
		
		if(args.length > 0){
			rC = Double.parseDouble(args[0]);
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		
		System.out.println("Lattice Energy Calculation of Alpha-phase Nitrogen");
		System.out.println("Using lattice sum method with truncation of " + rC + "A");
		System.out.println("with density of:" + density);
		
		HarmonicAlphaNitrogenModelLatticeSum test = new HarmonicAlphaNitrogenModelLatticeSum(Space3D.getInstance(3), numMolecule, density, rC);
		
	}
	
	private static final long serialVersionUID = 1L;
}
