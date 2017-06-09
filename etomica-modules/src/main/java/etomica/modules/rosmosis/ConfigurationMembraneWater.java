/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.IMoleculePositionDefinition;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;

public class ConfigurationMembraneWater implements Configuration {

    public ConfigurationMembraneWater(Simulation sim, Space _space) {
    	space = _space;
        soluteMoleFraction = 1;
        solutionChamberDensity = 0.5;
        solventChamberDensity = 0.5;
        numMembraneLayers = 2;
        membraneWidth = 4;
        this.sim = sim;
        positionDefinition = new MoleculePositionGeometricCenter(_space);
    }

    public void initializeCoordinates(Box box) {
        AtomActionTranslateBy translateBy = new AtomActionTranslateBy(space);
        Vector translationVector = translateBy.getTranslationVector();
        translationVector.E(0);
        MoleculeChildAtomAction translator = new MoleculeChildAtomAction(translateBy);
        
        box.setNMolecules(speciesSolute1, 0);
        box.setNMolecules(speciesSolute2, 0);
        box.setNMolecules(speciesSolvent, 0);
        box.setNMolecules(speciesMembrane, 0);
        
        Vector boxDimensions = box.getBoundary().getBoxSize();
        double boxLength = boxDimensions.getX(membraneDim);
        double chamberLength = 0.5 * boxLength - membraneTotalThickness;
        
        // solventChamber (middle, solvent-only)
        Box pretendBox = new Box(new BoundaryRectangularPeriodic(space, 1), space);
        sim.addBox(pretendBox);
        Vector pretendBoxDim = space.makeVector();
        pretendBoxDim.E(boxDimensions);
        pretendBoxDim.setX(membraneDim, chamberLength);
        pretendBox.getBoundary().setBoxSize(pretendBoxDim);
        int nMolecules = (int)Math.round(pretendBox.getBoundary().volume() * solventChamberDensity);
        System.out.println("adding "+nMolecules+" of species "+speciesSolvent);
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
        nSolutes = 2 * (nSolutes/2);
        pretendBox.setNMolecules(speciesSolute1, nSolutes);
//        pretendBox.setNMolecules(speciesSolute2, nSolutes/2);
        pretendBox.setNMolecules(speciesSolvent, nMolecules - nSolutes);
        configLattice.initializeCoordinates(pretendBox);
        // move molecules over to the real box
        ISpecies[] fluidSpecies = new ISpecies[]{speciesSolute1, speciesSolvent};
        for (int iSpecies=0; iSpecies<fluidSpecies.length; iSpecies++) {
            molecules = pretendBox.getMoleculeList(fluidSpecies[iSpecies]);
            for (int i=molecules.getMoleculeCount()-1; i>-1; i--) {
                // molecules will be reversed in order, but that's OK
                IMolecule molecule = molecules.getMolecule(i);
                pretendBox.removeMolecule(molecule);
                // we need to translate the molecules into the proper chamber
                double x = positionDefinition.position(molecule).getX(membraneDim);
                if (x < 0) {
                    translationVector.setX(membraneDim, -0.5*chamberLength - membraneTotalThickness);
                }
                else {
                    translationVector.setX(membraneDim, 0.5*chamberLength + membraneTotalThickness);
                }
                translator.actionPerformed(molecule);
                if (fluidSpecies[iSpecies] == speciesSolute1 && i % 2 == 0) {
                    // insert speciesSolute2 instead
                    IMolecule solute2 = speciesSolute2.makeMolecule();
                    translationVector.E(positionDefinition.position(molecule));
                    translator.actionPerformed(solute2);
                    molecule = solute2;
                    translationVector.E(0);
                }
                box.addMolecule(molecule);
            }
        }

        int pretendNumMembraneLayers = numMembraneLayers;
        double pretendMembraneThickness = membraneTotalThickness;
        double membraneThicknessPerLayer = membraneTotalThickness / numMembraneLayers;
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
        pretendBox.setNMolecules(speciesSolute1, 0);
        pretendBox.setNMolecules(speciesSolute2, 0);
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
                if (Math.abs(x - membraneCenter) > 0.5 * membraneTotalThickness) {
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

    /**
     * Returns the thickness (angstroms) of each layer of the membrane
     * (see getNumMembraneLayers)
     */
    public double getMembraneThickness() {
        return membraneTotalThickness;
    }

    public void setMembraneThickness(double newMembraneThickness) {
        membraneTotalThickness = newMembraneThickness;
    }

    /**
     * Returns the number of membrane layers.  If the membrane is in the yz
     * plane, then this is then number of layers in the x direction.
     * layers=2 means 1 unit cell.
     */
    public int getNumMembraneLayers() {
        return numMembraneLayers;
    }

    public void setNumMembraneLayers(int newNumMembraneLayers) {
        numMembraneLayers = newNumMembraneLayers;
    }

    /**
     * Returns the width of the membrane.  If the membrane is in the yz plane,
     * this is the number of lattice cells in the y direction and z direction.
     * membraneWidth = 2 means 2x2 cells.
     */
    public int getMembraneWidth() {
        return membraneWidth;
    }

    public void setMembraneWidth(int newMembraneWidth) {
        membraneWidth = newMembraneWidth;
    }

    /**
     * Returns the number density of the solvent chamber
     */
    public double getSolventChamberDensity() {
        return solventChamberDensity;
    }

    public void setSolventChamberDensity(double newSolventChamberDensity) {
        solventChamberDensity = newSolventChamberDensity;
    }

    /**
     * Returns the number density of the solution chamber
     */
    public double getSolutionChamberDensity() {
        return solutionChamberDensity;
    }

    public void setSolutionChamberDensity(double newSolutionChamberDensity) {
        solutionChamberDensity = newSolutionChamberDensity;
    }

    /**
     * Returns the mole fraction of solutes within the solute chamber.
     */
    public double getSoluteMoleFraction() {
        return soluteMoleFraction;
    }

    public void setSoluteMoleFraction(double newSoluteMoleFraction) {
        soluteMoleFraction = newSoluteMoleFraction;
    }

    /**
     * The dimension of the membrane (0=x, 1=y, 2=z).  x means that the
     * membrane exists in the yz plane.
     */
    public int getMembraneDim() {
        return membraneDim;
    }

    public void setMembraneDim(int newMembraneDim) {
        membraneDim = newMembraneDim;
    }

    public ISpecies getSpeciesSolute1() {
        return speciesSolute1;
    }

    public void setSpeciesSolute1(ISpecies newSpeciesSolute1) {
        speciesSolute1 = newSpeciesSolute1;
    }

    public ISpecies getSpeciesSolute2() {
        return speciesSolute2;
    }

    public void setSpeciesSolute2(ISpecies newSpeciesSolute2) {
        speciesSolute2 = newSpeciesSolute2;
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

    protected ISpecies speciesSolute1, speciesSolute2, speciesSolvent, speciesMembrane;
    protected double membraneTotalThickness;
    protected int numMembraneLayers, membraneWidth;
    protected double solventChamberDensity, solutionChamberDensity;
    protected double soluteMoleFraction;
    protected int membraneDim;
    protected final Simulation sim;
    private final Space space;
    protected IMoleculePositionDefinition positionDefinition;
}
