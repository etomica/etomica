package etomica.modules.rosmosis;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.simulation.ISimulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.species.Species;

public class ConfigurationMembraneWater implements Configuration {

    public ConfigurationMembraneWater(ISimulation sim) {
        soluteMoleFraction = 1;
        solutionChamberDensity = 0.5;
        solventChamberDensity = 0.5;
        numMembraneLayers = 2;
        membraneWidth = 4;
        this.sim = sim;
    }

    public void initializeCoordinates(Box box) {
        AtomActionTranslateBy translateBy = new AtomActionTranslateBy(sim.getSpace());
        IVector translationVector = translateBy.getTranslationVector();
        translationVector.E(0);
        AtomGroupAction translator = new AtomGroupAction(translateBy);
        
        box.setNMolecules(speciesSolute1, 0);
        box.setNMolecules(speciesSolute2, 0);
        box.setNMolecules(speciesSolvent, 0);
        box.setNMolecules(speciesMembrane, 0);
        
        IVector boxDimensions = box.getBoundary().getDimensions();
        double boxLength = boxDimensions.x(membraneDim);
        double membraneThickness = membraneThicknessPerLayer * numMembraneLayers;
        double chamberLength = 0.5 * boxLength - membraneThickness;
        
        // solventChamber (middle, solvent-only)
        Box pretendBox = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), null, 1), sim.getSpace());
        sim.addBox(pretendBox);
        IVector pretendBoxDim = sim.getSpace().makeVector();
        pretendBoxDim.E(boxDimensions);
        pretendBoxDim.setX(membraneDim, chamberLength);
        pretendBox.getBoundary().setDimensions(pretendBoxDim);
        int nMolecules = (int)Math.round(pretendBox.getBoundary().volume() * solventChamberDensity);
        System.out.println("adding "+nMolecules+" of species "+speciesSolvent);
        pretendBox.setNMolecules(speciesSolvent, nMolecules);
        ConfigurationLattice configLattice = new ConfigurationLattice(new LatticeCubicFcc(), sim.getSpace());
        configLattice.initializeCoordinates(pretendBox);
        // move molecules over to the real box
        AtomSet molecules = pretendBox.getMoleculeList(speciesSolvent);
        for (int i=nMolecules-1; i>-1; i--) {
            // molecules will be reversed in order, but that's OK
            IMolecule atom = (IMolecule)molecules.getAtom(i);
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
        Species[] fluidSpecies = new Species[]{speciesSolute1, speciesSolvent};
        for (int iSpecies=0; iSpecies<fluidSpecies.length; iSpecies++) {
            molecules = pretendBox.getMoleculeList(fluidSpecies[iSpecies]);
            for (int i=molecules.getAtomCount()-1; i>-1; i--) {
                // molecules will be reversed in order, but that's OK
                IMolecule atom = (IMolecule)molecules.getAtom(i);
                pretendBox.removeMolecule(atom);
                // we need to translate the molecules into the proper chamber
                double x = atom.getType().getPositionDefinition().position(atom).x(membraneDim);
                if (x < 0) {
                    translationVector.setX(membraneDim, -0.5*chamberLength - membraneThickness);
                }
                else {
                    translationVector.setX(membraneDim, 0.5*chamberLength + membraneThickness);
                }
                translator.actionPerformed(atom);
                if (fluidSpecies[iSpecies] == speciesSolute1 && i % 2 == 0) {
                    // insert speciesSolute2 instead
                    IMolecule solute2 = speciesSolute2.makeMolecule();
                    translationVector.E(atom.getType().getPositionDefinition().position(atom));
                    translator.actionPerformed(solute2);
                    atom = solute2;
                    translationVector.E(0);
                }
                box.addMolecule(atom);
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
        PrimitiveOrthorhombic primitive = new PrimitiveOrthorhombic(sim.getSpace());
        double a = boxDimensions.x(0) / membraneWidth;
        double b = boxDimensions.x(1) / membraneWidth;
        double c = boxDimensions.x(2) / membraneWidth;
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
        
        configLattice = new ConfigurationLattice(new BravaisLatticeCrystal(primitive, new BasisCubicFcc()), sim.getSpace());
        pretendBoxDim.E(boxDimensions);
        pretendBoxDim.setX(membraneDim, pretendMembraneThickness);
        pretendBox.getBoundary().setDimensions(pretendBoxDim);
        pretendBox.setNMolecules(speciesSolute1, 0);
        pretendBox.setNMolecules(speciesSolute2, 0);
        pretendBox.setNMolecules(speciesSolvent, 0);
        
        double[] shifts = new double[]{-0.25, 0.25};
        for (int iShift = 0; iShift<2; iShift++) {
            pretendBox.setNMolecules(speciesMembrane, nMolecules);
            configLattice.initializeCoordinates(pretendBox);
            
            double membraneShift = shifts[iShift]*boxDimensions.x(membraneDim) - membraneCenter;
            // move molecules over to the real box
            molecules = pretendBox.getMoleculeList(speciesMembrane);
            for (int i=molecules.getAtomCount()-1; i>-1; i--) {
                // molecules will be reversed in order, but that's OK
                IMolecule molecule = (IMolecule)molecules.getAtom(i);
                IAtomPositioned atom = (IAtomPositioned)molecule.getChildList().getAtom(0);
                double x = atom.getPosition().x(membraneDim);
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

    /**
     * Returns the thickness (angstroms) of each layer of the membrane
     * (see getNumMembraneLayers)
     */
    public double getMembraneThicknessPerLayer() {
        return membraneThicknessPerLayer;
    }

    public void setMembraneThicknessPerLayer(double newMembraneWidth) {
        membraneThicknessPerLayer = newMembraneWidth;
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

    public Species getSpeciesSolute1() {
        return speciesSolute1;
    }

    public void setSpeciesSolute1(Species newSpeciesSolute1) {
        speciesSolute1 = newSpeciesSolute1;
    }

    public Species getSpeciesSolute2() {
        return speciesSolute2;
    }

    public void setSpeciesSolute2(Species newSpeciesSolute2) {
        speciesSolute2 = newSpeciesSolute2;
    }

    public Species getSpeciesSolvent() {
        return speciesSolvent;
    }

    public void setSpeciesSolvent(Species newSpeciesSolvent) {
        speciesSolvent = newSpeciesSolvent;
    }

    public Species getSpeciesMembrane() {
        return speciesMembrane;
    }

    public void setSpeciesMembrane(Species newSpeciesMembrane) {
        speciesMembrane = newSpeciesMembrane;
    }

    protected Species speciesSolute1, speciesSolute2, speciesSolvent, speciesMembrane;
    protected double membraneThicknessPerLayer;
    protected int numMembraneLayers, membraneWidth;
    protected double solventChamberDensity, solutionChamberDensity;
    protected double soluteMoleFraction;
    protected int membraneDim;
    protected final ISimulation sim;
}
