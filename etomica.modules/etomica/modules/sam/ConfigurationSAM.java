package etomica.modules.sam;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.space.Space;

public class ConfigurationSAM implements Configuration {

    public ConfigurationSAM(ISimulation sim, Space space,
            ISpecies speciesMolecules, ISpecies speciesSurface) {
        this.sim = sim;
        this.space = space;
        this.speciesMolecules = speciesMolecules;
        this.speciesSurface = speciesSurface;
    }
    
    
    public void initializeCoordinates(IBox box) {
        Box pretendBox = new Box(box.getBoundary(), space);
        sim.addBox(pretendBox);
        box.setNMolecules(speciesMolecules, 0);
        int nMolecules = nCellsX*nCellsZ*basisMolecules.getScaledCoordinates().length;
        pretendBox.setNMolecules(speciesMolecules, nMolecules);
        
        IVector dim = box.getBoundary().getDimensions();
        dim.setX(0, nCellsX*cellSizeX);
        dim.setX(2, nCellsZ*cellSizeZ);
        double boxLengthY = dim.x(1);
        box.setDimensions(dim);
        pretendBox.setDimensions(dim);
        Primitive primitive = new PrimitiveOrthorhombic(space, cellSizeX, boxLengthY, cellSizeZ);
        BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basisMolecules);
        Configuration config = new ConfigurationLattice(lattice, space);
        config.initializeCoordinates(pretendBox);

        IAtomSet molecules = pretendBox.getMoleculeList(speciesMolecules);
        double y0 = ((IAtomPositioned)((IMolecule)molecules.getAtom(0)).getChildList().getAtom(0)).getPosition().x(1);
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = (IMolecule)molecules.getAtom(0);
            pretendBox.removeMolecule(molecule);
            box.addMolecule(molecule);
        }

        nMolecules = nCellsX*nCellsZ*basisSurface.getScaledCoordinates().length;
        pretendBox.setNMolecules(speciesMolecules, 0);
        pretendBox.setNMolecules(speciesSurface, nMolecules);
        lattice = new BravaisLatticeCrystal(primitive, basisSurface);
        config = new ConfigurationLattice(lattice, space);
        config.initializeCoordinates(pretendBox);

        molecules = pretendBox.getMoleculeList(speciesSurface);
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = (IMolecule)molecules.getAtom(0);
            pretendBox.removeMolecule(molecule);
            box.addMolecule(molecule);
            IVector pos = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
            pos.setX(1, y0-yOffset);
        }
    }

    public double getYOffset() {
        return yOffset;
    }
    
    public void setYOffset(double newYOffset) {
        yOffset = newYOffset;
    }
    
    public double getCellSizeX() {
        return cellSizeX;
    }


    public void setCellSizeX(double cellSizeX) {
        this.cellSizeX = cellSizeX;
    }


    public double getCellSizeZ() {
        return cellSizeZ;
    }


    public void setCellSizeZ(double cellSizeZ) {
        this.cellSizeZ = cellSizeZ;
    }


    public int getNCellsX() {
        return nCellsX;
    }


    public void setNCellsX(int cellsX) {
        nCellsX = cellsX;
    }


    public int getNCellsZ() {
        return nCellsZ;
    }


    public void setNCellsZ(int cellsZ) {
        nCellsZ = cellsZ;
    }


    public Basis getBasisMolecules() {
        return basisMolecules;
    }


    public void setBasisMolecules(Basis basisMolecules) {
        this.basisMolecules = basisMolecules;
    }


    public Basis getBasisSurface() {
        return basisSurface;
    }


    public void setBasisSurface(Basis basisSurface) {
        this.basisSurface = basisSurface;
    }
    
    protected final Space space;
    protected final ISimulation sim;
    protected final ISpecies speciesMolecules;
    protected final ISpecies speciesSurface;
    protected double cellSizeX, cellSizeZ;
    protected int nCellsX, nCellsZ;
    protected Basis basisMolecules;
    protected Basis basisSurface;
    protected double yOffset;
}
