package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * Class to create a scaled basis for coordinate definition
 * 	
 * 
 * @author Tai Boon Tan
 *
 */
public class BasisNMOneCell {
	
	public BasisNMOneCell(ISpace space, int numAtoms, double density){
		this.space = space;
		this.numAtoms = numAtoms;
		this.density = density;
	}
	
	public IVectorMutable[] getScaledBasis(){
		
		Simulation simulation = new Simulation(space);
        SpeciesSpheresMono speciesImag = new SpeciesSpheresMono(simulation, space);
        simulation.addSpecies(speciesImag);
        
        Box boxImag = new Box(space);
        simulation.addBox(boxImag);
        boxImag.setNMolecules(speciesImag, numAtoms);
        
        double L = Math.pow(4.0/density, 1.0/3.0);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        Primitive primitiveImag = new PrimitiveCubic(space, L);
        Boundary boundaryImag = new BoundaryRectangularPeriodic(space, n*L);
        Basis basisImag = new BasisCubicFcc();
        boxImag.setBoundary(boundaryImag);
        
        BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitiveImag, basisImag);
                
        ConfigurationLattice configLattice = new ConfigurationLattice(lattice, space);
        configLattice.initializeCoordinates(boxImag);
        IVectorMutable halfBoxLength = (IVectorMutable)boxImag.getBoundary().getBoxSize();
        halfBoxLength.TE(0.5);
        
        IVectorMutable[] basisVector = new IVectorMutable[numAtoms];
              
        
        for (int i=0; i<numAtoms; i++){
        	basisVector[i] = space.makeVector();
        	IAtom atom = boxImag.getLeafList().getAtom(i);
        	
        	IVectorMutable pos = atom.getPosition();
        	pos.PE(halfBoxLength);
        	pos.TE(1/(n*L));
        	basisVector[i].E(pos);
        	
        }
		
		return basisVector;
		
	}
	
	protected ISpace space;
	protected int numAtoms;
	protected double density;
	
}
