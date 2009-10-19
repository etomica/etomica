package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.config.ConfigurationLatticeSimple;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.species.SpeciesSpheresMono;

/**
 * Basis for a cell constructed using multiple cells.
 * 
 * @author Andrew Schultz
 */
public class BasisBigCell extends Basis {

    private static final long serialVersionUID = 1L;

    public BasisBigCell(ISpace space, Primitive primitive, Basis subBasis, int[] nSubCells) {
        super(makeScaledCoordinates(space, primitive, subBasis, nSubCells));
    }

    protected static IVector[] makeScaledCoordinates(ISpace space, Primitive primitive, Basis subBasis, int[] nSubCells) {
        // make pretend sim, species and box so we can find the appropriate coordinates
        ISimulation sim = new Simulation(space);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.getSpeciesManager().addSpecies(species);
        IBoundary boundary = new BoundaryDeformableLattice(primitive, nSubCells);
        IBox box = new Box(boundary, space);
        sim.addBox(box);
        int numMolecules = subBasis.getScaledCoordinates().length;
        for (int i=0; i<nSubCells.length; i++) {
            numMolecules *= nSubCells[i];
        }
        box.setNMolecules(species, numMolecules);
        ConfigurationLatticeSimple configLattice = new ConfigurationLatticeSimple(new BravaisLatticeCrystal(primitive, subBasis), space);
        configLattice.initializeCoordinates(box);

        // determine the transformation from real coordinates to scaled coordinates
        IVectorMutable[] edges = new IVectorMutable[space.D()];
        for (int i=0; i<edges.length; i++) {
            edges[i] = space.makeVector();
            edges[i].E(boundary.getEdgeVector(i));
        }
        Tensor boundaryTensor = space.makeTensor();
        boundaryTensor.E(edges);
        boundaryTensor.invert();

        // retrieve real coordinates and scale them
        IAtomList atomList = box.getLeafList();
        IVectorMutable[] pos = new IVectorMutable[atomList.getAtomCount()];
        for (int i=0; i<atomList.getAtomCount(); i++) {
            pos[i] = space.makeVector();
            pos[i].E(atomList.getAtom(i).getPosition());
            boundaryTensor.transform(pos[i]);
            // coordinates now range from -0.5 to +0.5, we want 0.0 to 1.0
            pos[i].PE(0.5);
        }
        return pos;
    }
}
