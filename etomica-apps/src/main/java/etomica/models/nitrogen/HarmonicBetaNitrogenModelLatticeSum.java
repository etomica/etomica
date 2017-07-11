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
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.models.nitrogen.LatticeSumCrystalMolecular.DataGroupLSC;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.BasisBigCell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Energy;
import etomica.units.Joule;


/**
 *  Lattice sum class for Beta-phase Nitrogen
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class HarmonicBetaNitrogenModelLatticeSum extends Simulation{

	
	public HarmonicBetaNitrogenModelLatticeSum(Space space, int numMolecule, double density, double rC) {
		super(space);
						
	  	double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		//System.out.println("aDim: " + aDim + " ;cDim: " + cDim);
		
		int [] nCells = new int[]{1,2,1};
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, nCells);
        
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		SpeciesN2 species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		SpeciesN2B ghostSpecies = new SpeciesN2B(space);
		ghostSpecies.setConformation(conformation);
		addSpecies(ghostSpecies);
		
		Box box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		Box ghostBox = new Box(space);
		addBox(ghostBox);
		ghostBox.setNMolecules(ghostSpecies, 1);
		
		Primitive primitive = new PrimitiveTriclinic(space, aDim, 2*aDim, cDim, Math.PI*(90/180.0),Math.PI*(90/180.0),Math.PI*(120/180.0));

    	BetaPhaseLatticeParameterLS parameters = new BetaPhaseLatticeParameterLS();
		double[][] param = parameters.getParameter(density);
		
		CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBetaLatticeSum();
		coordinateDef.setIsDoLatticeSum();
		coordinateDef.setOrientationVectorBetaLatticeSum(space, density, param);
		coordinateDef.initializeCoordinates(new int[]{1,1,1});
	
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
		
		int nLayer = (int)Math.round(rC/aDim + 0.5);
		
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
        System.out.println(density+ " lattice energy [sim unit]:  " + latEnergy + " ;[kJ/mol]: " + Joule.UNIT.fromSim(latEnergy)*avogradoConst/1000);

	}
	
	public static void main (String[] args){
		
		int numMolecule =4;
		double density = 0.0230;
		double rC = 80;
	
		if(args.length > 0){
			rC = Double.parseDouble(args[0]);
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		System.out.println("Lattice Energy Calculation of Beta-phase Nitrogen");
		System.out.println("Using lattice sum method with truncation of " + rC + "A");
		System.out.println("with density of:" + density);
		
		HarmonicBetaNitrogenModelLatticeSum sim = new HarmonicBetaNitrogenModelLatticeSum(Space3D.getInstance(3), numMolecule, density, rC);
		
	}
	
	private static final long serialVersionUID = 1L;
}
