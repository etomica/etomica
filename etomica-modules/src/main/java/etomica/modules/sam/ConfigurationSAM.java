/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLatticeSimple;
import etomica.config.ConformationChainZigZag;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.PotentialMasterList;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

public class ConfigurationSAM implements Configuration {

    public ConfigurationSAM(Simulation sim, Space space,
                            ISpecies speciesMolecules, ISpecies speciesSurface,
                            PotentialMasterList potentialMaster) {
        this.sim = sim;
        this.space = space;
        this.speciesMolecules = speciesMolecules;
        this.speciesSurface = speciesSurface;
        moleculeOffset = space.makeVector();
        this.potentialMaster = potentialMaster;
        conformation = new ConformationChainZigZag[4];
    }
    
    public ConformationChainZigZag getConformation(int iChain) {
        return conformation[iChain];
    }
    
    public void setConformation(int iChain, ConformationChainZigZag newConformation) {
        conformation[iChain] = newConformation;
    }
    
    public void setMoleculeOffset(Vector newMoleculeOffset) {
        moleculeOffset.E(newMoleculeOffset);
    }
    
    public Vector getMoleculeOffset() {
        return moleculeOffset;
    }
    
    public void initializeCoordinates(Box box) {
        Box pretendBox = new Box(box.getBoundary(), space);
        sim.addBox(pretendBox);
        potentialMaster.getNbrCellManager(pretendBox).setDoApplyPBC(true);
        potentialMaster.getNeighborManager(pretendBox).setDoApplyPBC(false);
        box.setNMolecules(speciesMolecules, 0);
        box.setNMolecules(speciesSurface, 0);
        int nMolecules = nCellsX*nCellsZ*basisMolecules.getScaledCoordinates().length;
        pretendBox.setNMolecules(speciesMolecules, nMolecules);
        
        Vector dim = space.makeVector();
        dim.E(box.getBoundary().getBoxSize());
        dim.setX(0, nCellsX*cellSizeX);
        dim.setX(2, nCellsZ*cellSizeZ);
        double boxLengthY = dim.getX(1);
        box.getBoundary().setBoxSize(dim);
        pretendBox.getBoundary().setBoxSize(dim);
        Primitive primitive = new PrimitiveOrthorhombic(space, cellSizeX, boxLengthY, cellSizeZ);
        BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basisMolecules);
        ConfigurationLatticeSimple config = new ConfigurationLatticeSimple(lattice, space);
        config.initializeCoordinates(pretendBox);
        
        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
        translator.getTranslationVector().E(moleculeOffset);
        MoleculeChildAtomAction groupTranslator = new MoleculeChildAtomAction(translator);
        
        Vector offset = space.makeVector();

        IMoleculeList molecules = pretendBox.getMoleculeList(speciesMolecules);
        double y0 = molecules.get(0).getChildList().get(0).getPosition().getX(1) + moleculeOffset.getX(1);
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.get(0);
            pretendBox.removeMolecule(molecule);
            groupTranslator.actionPerformed(molecule);
            box.addMolecule(molecule);
            if (i % 2 == 0 && conformation[1] != null) {
                if ((i-1)/(nCellsZ*2) % 2 == 1 && conformation[3] != null) {
                    offset.E(molecule.getChildList().get(0).getPosition());
                    conformation[3].initializePositions(molecule.getChildList());
                    offset.ME(molecule.getChildList().get(0).getPosition());
                    translator.getTranslationVector().E(offset);
                    groupTranslator.actionPerformed(molecule);
                    translator.getTranslationVector().E(moleculeOffset);        
                }
                else {
                    offset.E(molecule.getChildList().get(0).getPosition());
                    conformation[1].initializePositions(molecule.getChildList());
                    offset.ME(molecule.getChildList().get(0).getPosition());
                    translator.getTranslationVector().E(offset);
                    groupTranslator.actionPerformed(molecule);
                    translator.getTranslationVector().E(moleculeOffset);
                }
            }
            else if ((i-1)/(nCellsZ*2) % 2 == 1 && conformation[2] != null) {
                offset.E(molecule.getChildList().get(0).getPosition());
                conformation[2].initializePositions(molecule.getChildList());
                offset.ME(molecule.getChildList().get(0).getPosition());
                translator.getTranslationVector().E(offset);
                groupTranslator.actionPerformed(molecule);
                translator.getTranslationVector().E(moleculeOffset);
            }
        }

        nMolecules = nCellsX*nCellsZ*basisSurface.getScaledCoordinates().length;
        pretendBox.setNMolecules(speciesMolecules, 0);
        pretendBox.setNMolecules(speciesSurface, nMolecules);
        lattice = new BravaisLatticeCrystal(primitive, basisSurface);
        config = new ConfigurationLatticeSimple(lattice, space);
        config.initializeCoordinates(pretendBox);

        molecules = pretendBox.getMoleculeList(speciesSurface);
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.get(0);
            pretendBox.removeMolecule(molecule);
            box.addMolecule(molecule);
            Vector pos = molecule.getChildList().get(0).getPosition();
            pos.setX(1, y0-yOffset);
        }
        sim.removeBox(pretendBox);
    }

    public double getSurfaceYOffset() {
        return yOffset;
    }
    
    public void setSurfaceYOffset(double newYOffset) {
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
    protected final Simulation sim;
    protected final ISpecies speciesMolecules;
    protected final ISpecies speciesSurface;
    protected double cellSizeX, cellSizeZ;
    protected int nCellsX, nCellsZ;
    protected Basis basisMolecules;
    protected Basis basisSurface;
    protected double yOffset;
    protected final Vector moleculeOffset;
    protected ConformationChainZigZag[] conformation;
    protected PotentialMasterList potentialMaster;
}
