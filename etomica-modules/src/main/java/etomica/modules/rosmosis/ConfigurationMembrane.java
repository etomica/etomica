/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;


import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.atom.IMoleculePositionDefinition;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;

public class ConfigurationMembrane implements Configuration {

    public ConfigurationMembrane(Simulation sim, Space _space) {
        soluteMoleFraction = 1;
        solutionChamberDensity = 0.5;
        solventChamberDensity = 0.5;
        numMembraneLayers = 2;
        membraneWidth = 4;
        this.sim = sim;
        this.space = _space;
        positionDefinition = new MoleculePositionGeometricCenter(space);
    }

    public void initializeCoordinates(Box box) {
        AtomActionTranslateBy translateBy = new AtomActionTranslateBy(space);
        Vector translationVector = translateBy.getTranslationVector();
        MoleculeChildAtomAction translator = new MoleculeChildAtomAction(translateBy);
        translationVector.E(0);
        
        box.setNMolecules(speciesSolute, 0);
        box.setNMolecules(speciesSolvent, 0);
        box.setNMolecules(speciesMembrane, 0);
        
        Vector boxDimensions = box.getBoundary().getBoxSize();
        double boxLength = boxDimensions.getX(membraneDim);
        double membraneThickness = membraneThicknessPerLayer * numMembraneLayers;
        double chamberLength = 0.5 * boxLength - membraneThickness;
        
        // solventChamber (middle, solvent-only)
        Box pretendBox = new Box(new BoundaryRectangularPeriodic(space, 1), space);
        sim.addBox(pretendBox);
        Vector pretendBoxDim = space.makeVector();
        pretendBoxDim.E(boxDimensions);
        pretendBoxDim.setX(membraneDim, chamberLength);
        pretendBox.getBoundary().setBoxSize(pretendBoxDim);
        int nMolecules = (int)Math.round(pretendBox.getBoundary().volume() * solventChamberDensity);
        pretendBox.setNMolecules(speciesSolvent, nMolecules);
        ConfigurationLattice configLattice = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configLattice.initializeCoordinates(pretendBox);
        // move molecules over to the real box
        IMoleculeList molecules = pretendBox.getMoleculeList(speciesSolvent);
        for (int i=nMolecules-1; i>-1; i--) {
            // molecules will be reversed in order, but that's OK
            IMolecule atom = molecules.getMolecule(i);
            pretendBox.removeMolecule(atom);
            box.addMolecule(atom);
        }

        nMolecules = (int)Math.round(pretendBox.getBoundary().volume() * solutionChamberDensity);
        int nSolutes = (int)(nMolecules * soluteMoleFraction);
        pretendBox.setNMolecules(speciesSolute, nSolutes);
        pretendBox.setNMolecules(speciesSolvent, nMolecules - nSolutes);
        configLattice.initializeCoordinates(pretendBox);
        // move molecules over to the real box
        ISpecies[] fluidSpecies = new ISpecies[]{speciesSolute, speciesSolvent};
        for (int iSpecies=0; iSpecies<2; iSpecies++) {
            molecules = pretendBox.getMoleculeList(fluidSpecies[iSpecies]);
            for (int i=molecules.getMoleculeCount()-1; i>-1; i--) {
                // molecules will be reversed in order, but that's OK
                IMolecule molecule = molecules.getMolecule(i);
                pretendBox.removeMolecule(molecule);
                // we need to translate the molecules into the proper chamber
                double x = positionDefinition.position(molecule).getX(membraneDim);
                if (x < 0) {
                    translationVector.setX(membraneDim, -0.5*chamberLength - membraneThickness);
                }
                else {
                    translationVector.setX(membraneDim, 0.5*chamberLength + membraneThickness);
                }
                translator.actionPerformed(molecule);
                box.addMolecule(molecule);
            }
        }
        
        int pretendNumMembraneLayers = numMembraneLayers;
        double pretendMembraneThickness = membraneThickness;
        double membraneCenter = 0;
        if (numMembraneLayers % 2 == 1) {
            // we want an odd number of layers, which ConfigurationLattice can't
            // handle.  So add a pretend layer of atoms, which we'll drop later.
            pretendMembraneThickness /= numMembraneLayers;
            pretendNumMembraneLayers++;
            pretendMembraneThickness *= pretendNumMembraneLayers;
            membraneCenter = -0.5 * membraneThicknessPerLayer;
        }
        
        nMolecules = 2*membraneWidth * membraneWidth * pretendNumMembraneLayers;
        PrimitiveOrthorhombic primitive = new PrimitiveOrthorhombic(space);
        double a = boxDimensions.getX(0) / membraneWidth;
        double b = boxDimensions.getX(1) / membraneWidth;
        double c = boxDimensions.getX(2) / membraneWidth;
        switch (membraneDim) {
            case 0:
                a = 2*pretendMembraneThickness / pretendNumMembraneLayers;
                break;
            case 1:
                b = 2*pretendMembraneThickness / pretendNumMembraneLayers;
                break;
            case 2:
                c = 2*pretendMembraneThickness / pretendNumMembraneLayers;
                break;
        }
        primitive.setSizeA(a);
        primitive.setSizeB(b);
        primitive.setSizeC(c);
        
        configLattice = new ConfigurationLattice(new BravaisLatticeCrystal(primitive, new BasisCubicFcc()), space);
        pretendBoxDim.E(boxDimensions);
        pretendBoxDim.setX(membraneDim, pretendMembraneThickness);
        pretendBox.getBoundary().setBoxSize(pretendBoxDim);
        pretendBox.setNMolecules(speciesSolute, 0);
        pretendBox.setNMolecules(speciesSolvent, 0);
        
        double[] shifts = new double[]{-0.25, 0.25};
        for (int iShift = 0; iShift<2; iShift++) {
            pretendBox.setNMolecules(speciesMembrane, nMolecules);
            configLattice.initializeCoordinates(pretendBox);
            
            double membraneShift = shifts[iShift]*boxDimensions.getX(membraneDim) - membraneCenter;
            // move molecules over to the real box
            molecules = pretendBox.getMoleculeList(speciesMembrane);
            for (int i=molecules.getMoleculeCount()-1; i>-1; i--) {
                // molecules will be reversed in order, but that's OK
                IMolecule molecule = molecules.getMolecule(i);
                IAtom atom = molecule.getChildList().getAtom(0);
                double x = atom.getPosition().getX(membraneDim);
                if (Math.abs(x - membraneCenter) > 0.5 * membraneThickness) {
                    // we encountered a pretend atom in our pretend box!
                    continue;
                }
                atom.getPosition().setX(membraneDim, x + membraneShift);
                pretendBox.removeMolecule(molecule);
                box.addMolecule(molecule);
            }
        }
        
        sim.removeBox(pretendBox);
    }

    public double getMembraneThicknessPerLayer() {
        return membraneThicknessPerLayer;
    }

    public void setMembraneThicknessPerLayer(double newMembraneWidth) {
        membraneThicknessPerLayer = newMembraneWidth;
    }

    public int getNumMembraneLayers() {
        return numMembraneLayers;
    }

    public void setNumMembraneLayers(int newNumMembraneLayers) {
        numMembraneLayers = newNumMembraneLayers;
    }

    public int getMembraneWidth() {
        return membraneWidth;
    }

    public void setMembraneWidth(int newMembraneWidth) {
        membraneWidth = newMembraneWidth;
    }

    public double getSolventChamberDensity() {
        return solventChamberDensity;
    }

    public void setSolventChamberDensity(double newSolventChamberDensity) {
        solventChamberDensity = newSolventChamberDensity;
    }

    public double getSolutionChamberDensity() {
        return solutionChamberDensity;
    }

    public void setSolutionChamberDensity(double newSolutionChamberDensity) {
        solutionChamberDensity = newSolutionChamberDensity;
    }

    public double getSoluteMoleFraction() {
        return soluteMoleFraction;
    }

    public void setSoluteMoleFraction(double newSoluteMoleFraction) {
        soluteMoleFraction = newSoluteMoleFraction;
    }

    public int getMembraneDim() {
        return membraneDim;
    }

    public void setMembraneDim(int newMembraneDim) {
        membraneDim = newMembraneDim;
    }

    public ISpecies getSpeciesSolute() {
        return speciesSolute;
    }

    public void setSpeciesSolute(ISpecies newSpeciesSolute) {
        speciesSolute = newSpeciesSolute;
    }

    public ISpecies getSpeciesSolvent() {
        return speciesSolvent;
    }

    public void setSpeciesSolvent(ISpecies newSpeciesSolvent) {
        speciesSolvent = newSpeciesSolvent;
    }

    public ISpecies getSpeciesMembrane() {
        return speciesMembrane;
    }

    public void setSpeciesMembrane(ISpecies newSpeciesMembrane) {
        speciesMembrane = newSpeciesMembrane;
    }

    protected ISpecies speciesSolute, speciesSolvent, speciesMembrane;
    protected double membraneThicknessPerLayer;
    protected int numMembraneLayers, membraneWidth;
    protected double solventChamberDensity, solutionChamberDensity;
    protected double soluteMoleFraction;
    protected int membraneDim;
    protected final Simulation sim;
    private final Space space;
    protected IMoleculePositionDefinition positionDefinition;
}
