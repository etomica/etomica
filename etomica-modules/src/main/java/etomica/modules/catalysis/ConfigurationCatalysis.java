/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.atom.AtomLeafAgentManager;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationLatticeSimple;
import etomica.config.ConformationChainZigZag;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;

/**
 * Configuration for catalysis module.  Places molecules within the box on an
 * FCC lattice, and the surface atoms at the bottom of the box in a hexagonal
 * lattice.
 * 
 * @author Andrew Schultz
 */
public class ConfigurationCatalysis implements Configuration {

    public ConfigurationCatalysis(Simulation sim, Space space,
                                  ISpecies speciesSurface, ISpecies speciesC, ISpecies speciesO, AtomLeafAgentManager agentManager) {
        this.sim = sim;
        this.space = space;
        this.speciesSurface = speciesSurface;
        this.speciesO = speciesO;
        this.speciesC = speciesC;
        this.agentManager = agentManager;
        moleculeOffset = space.makeVector();
        conformation = new ConformationChainZigZag[4];
    }
    
    public void setMoleculeOffset(Vector newMoleculeOffset) {
        moleculeOffset.E(newMoleculeOffset);
    }
    
    public Vector getMoleculeOffset() {
        return moleculeOffset;
    }
    
    public void initializeCoordinates(Box box) {
        box.setNMolecules(speciesSurface, 0);
        box.setNMolecules(speciesO, 0);
        box.setNMolecules(speciesC, 0);
        
        Vector dim = space.makeVector();
        dim.E(box.getBoundary().getBoxSize());
        dim.setX(0, nCellsX*cellSizeX);
        dim.setX(1, 0.9*dim.getX(1));
        dim.setX(2, nCellsZ*cellSizeZ);
        box.getBoundary().setBoxSize(dim);
        
        Box pretendBox = new Box(box.getBoundary(), space);
        sim.addBox(pretendBox);
        
        //initialize the "molecules"
        pretendBox.setNMolecules(speciesO, nO2);
        pretendBox.setNMolecules(speciesC, nCO);
        LatticeCubicFcc lattice = new LatticeCubicFcc(space);
        Configuration config = new ConfigurationLattice(lattice, space);
        config.initializeCoordinates(pretendBox);
        dim.setX(1, dim.getX(1)/0.9);
        box.getBoundary().setBoxSize(dim);
        
        IMoleculeList molecules = pretendBox.getMoleculeList();
        Vector shift = space.makeVector();
        shift.setX(0, -1.901);
        while (molecules.getMoleculeCount()>0) {
            IMolecule molecule1 = molecules.getMolecule(0);
            pretendBox.removeMolecule(molecule1);
            box.addMolecule(molecule1);
            IMolecule molecule2 = speciesO.makeMolecule();
            box.addMolecule(molecule2);
            IAtom atom1 = molecule1.getChildList().getAtom(0);
            Vector pos1 = atom1.getPosition();
            IAtom atom2 = molecule2.getChildList().getAtom(0);
            Vector pos2 = atom2.getPosition();
            pos2.Ev1Mv2(pos1, shift);
            pos1.PE(shift);
            ((CatalysisAgent)agentManager.getAgent(atom1)).bondedAtom1 = atom2;
            ((CatalysisAgent)agentManager.getAgent(atom2)).bondedAtom1 = atom1;
        }

        Primitive primitive = new PrimitiveOrthorhombic(space, cellSizeX, dim.getX(1), cellSizeZ);
        BasisOrthorhombicHexagonal3D basisSurface = new BasisOrthorhombicHexagonal3D(space);
        int nMolecules = nCellsX*nCellsZ*basisSurface.getScaledCoordinates().length;
        pretendBox.setNMolecules(speciesSurface, nMolecules);
        BravaisLatticeCrystal latticeSurface = new BravaisLatticeCrystal(primitive, basisSurface);
        config = new ConfigurationLatticeSimple(latticeSurface, space);
        config.initializeCoordinates(pretendBox);

        molecules = pretendBox.getMoleculeList(speciesSurface);
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.getMolecule(0);
            pretendBox.removeMolecule(molecule);
            box.addMolecule(molecule);
        }
        
//        pretendBox.setNMolecules(speciesSurface, nMolecules);
//        config.initializeCoordinates(pretendBox);
//
//        molecules = pretendBox.getMoleculeList(speciesSurface);
//        for (int i=0; i<nMolecules; i++) {
//            IMolecule molecule = molecules.getMolecule(0);
//            pretendBox.removeMolecule(molecule);
//            IAtom atom = molecule.getChildList().getAtom(0);
//            IVectorMutable pos = ((IAtomPositioned)atom).getPosition();
//            pos.setX(1, -pos.getX(1)-0.1);
//            box.addMolecule(molecule);
//        }

        sim.removeBox(pretendBox);
    }
    
    public void setNumCO(int newNumCO) {
        nCO = newNumCO;
    }
    
    public void setNumO2(int newNumO2) {
        nO2 = newNumO2;
    }
    
    public int getNumCO() {
        return nCO;
    }

    public int getNumO2() {
        return nO2;
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

    protected final Space space;
    protected final Simulation sim;
    protected final ISpecies speciesSurface, speciesO, speciesC;
    protected int nCO, nO2;
    protected final AtomLeafAgentManager agentManager;
    protected double cellSizeX, cellSizeZ;
    protected int nCellsX, nCellsZ;
    protected double yOffset;
    protected final Vector moleculeOffset;
    protected ConformationChainZigZag[] conformation;
    protected PotentialMasterList potentialMaster;
    
    
    public static class BasisOrthorhombicHexagonal3D extends Basis {
        /**
         * Makes a fcc 4-atom basis.
         */
        public BasisOrthorhombicHexagonal3D(Space space) {
            super(makeScaledPositions(space));
        }
        
        private static final Vector[] makeScaledPositions(Space space) {
            Vector[] v = new Vector[2];
            for (int i=0; i<2; i++) {
                v[i] = space.makeVector();
            }
            v[0].E(0);
            v[1].setX(0, 0.5);
            v[1].setX(2, 0.5);
            return v;
        };
        
        private static final long serialVersionUID = 1L;
    }

}
