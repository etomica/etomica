package etomica.box;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.vecarray.VecArray;

import java.util.AbstractList;
import java.util.List;
import java.util.RandomAccess;

public class BoxFast {
    private final Boundary boundary;
    private final List<ISpecies> speciesList;

    private int atomCount;
    private int moleculeCount;

    private final VecArray atomPositions;
    private final VecArray atomVelocities;

    // For every species, an array that has the atom indices for each molecule.
    // Since number of atoms in a species is known, the indices are packed in one array
    private int[][] moleculeAtoms;

    public BoxFast(List<ISpecies> speciesList) {
        this(new BoundaryRectangularPeriodic(Space3D.getInstance()), speciesList);
    }

    public BoxFast(Boundary boundary, List<ISpecies> speciesList) {
        this.boundary = boundary;
        this.speciesList = speciesList;
    }

    public Space getSpace() { return Space3D.getInstance(); }

    public void setNMolecules(ISpecies species, int n) {

    }

    private void addMolecule(IMolecule molecule) {
        
    }

    private static class FakeAtom implements IAtomKinetic {

        @Override
        public Vector getVelocity() {
            return null;
        }

        @Override
        public int getIndex() {
            return 0;
        }

        @Override
        public void setIndex(int index) {

        }

        @Override
        public int getLeafIndex() {
            return 0;
        }

        @Override
        public void setLeafIndex(int newGlobalIndex) {

        }

        @Override
        public void setParent(IMolecule newParent) {

        }

        @Override
        public IMolecule getParentGroup() {
            return null;
        }

        @Override
        public AtomType getType() {
            return null;
        }

        @Override
        public Vector getPosition() {
            return null;
        }
    }

    private static class FakeAtomList extends AbstractList<IAtom> implements IAtomList, RandomAccess {

        @Override
        public IAtom get(int index) {
            return null;
        }

        @Override
        public int size() {
            return 0;
        }
    }
}
