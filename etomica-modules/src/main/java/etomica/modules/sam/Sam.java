/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationChainZigZag;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.IOrientation3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Calorie;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.Constants;

/**
 * Self-assembled monolayer module.
 * @author Andrew Schultz
 *
 */
public class Sam extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final P1Sinusoidal p1SurfaceBond;
    public final double sinusoidalB;
    public final P2Harmonic p2SurfaceBond;
    public final double harmonicStrength;
    public SpeciesAlkaneThiol species;
    public SpeciesSpheresMono speciesSurface;
    public Box box;
    public IntegratorVelocityVerletSAM integrator;
    public ActivityIntegrate activityIntegrate;
    public IMoleculePositionDefinition positionDefinition;
    public P1WCAWall wallPotential;
    public ConfigurationSAM config;
    public P2LennardJones p2CH2, p2CH3, p2CH2CH3, p2S, p2SCH2;
    public P2SoftSphericalTruncatedSwitched p2CH2t, p2CH3t, p2CH2CH3t, p2St, p2SCH2t;
    public P2Harmonic p2BondCC, p2BondCS;
    public P3BondAngle p3Bond;
    public P4BondTorsion p4BondCCCC, p4BondCCCS;
    public PotentialGroup p1Intra;
    public P2SoftSphericalTruncatedSwitched p2SulfurSurfaceLJ, p2CH2Surface;
    public CriterionTether3 criterion3;
    public int chainLength;
    public double chainTheta, chainPsi;
    public double[] chainPhi;
    public double bondL_CC;
    public double bondTheta;
    public PotentialCalculationForceSumWall forceSum;
    public PotentialMasterList potentialMaster;
    public int numZCells, numXCells;
    public double sizeCellZ, sizeCellX;
    public double sigmaCH2;

    public Sam() {
        super(Space.getInstance(3));

        //species and potentials
        chainLength = 16;
        species = new SpeciesAlkaneThiol(space, chainLength - 1);
        species.setIsDynamic(true);
        addSpecies(species);

        speciesSurface = new SpeciesSpheresMono(this, space);
        speciesSurface.setIsDynamic(true);
        ((ElementSimple) speciesSurface.getLeafType().getElement()).setMass(Double.POSITIVE_INFINITY);
        addSpecies(speciesSurface);

        positionDefinition = new MoleculePositionGeometricCenter(space);
        sigmaCH2 = 3.95;
        potentialMaster = new PotentialMasterList(this, 2.8 * sigmaCH2, space); //List(this, 2.0);

        numXCells = 4;
        numZCells = 2;
        // gold has FCC unit cell, a=4.0782A
        sizeCellZ = 4.0782 / Math.sqrt(2) * 3; //Math.sqrt(3)*sizeCellX;
        sizeCellX = sizeCellZ / Math.sqrt(3);

        //construct box
        box = this.makeBox(new BoundaryRectangularSlit(1, space));
        Vector dim = space.makeVector();
        dim.E(new double[]{sizeCellX * numXCells, chainLength * 2.4, sizeCellZ * numZCells});
        box.getBoundary().setBoxSize(dim);

        bondL_CC = 1.54;
        double bondL_CS = 1.82;
        bondTheta = Math.PI * 114 / 180;
        ConformationChainZigZag conformation = new ConformationChainZigZag(space);
        species.setConformation(conformation);
        chainTheta = 0;
        chainPsi = 0;
        chainPhi = new double[4];

        config = new ConfigurationSAM(this, space, species, speciesSurface, potentialMaster);
        Basis alkaneBasis = new Basis(new Vector[]{space.makeVector(new double[]{1.0 / 6.0, 0, 1.0 / 6.0}), space.makeVector(new double[]{2.0 / 3.0, 0, 2.0 / 3.0})});
        Basis surfaceBasis = new Basis(new Vector[]{
                space.makeVector(new double[]{2.0 / 6.0, 0, 0}),
                space.makeVector(new double[]{5.0 / 6.0, 0, 1.0 / 6.0}),
                space.makeVector(new double[]{2.0 / 6.0, 0, 2.0 / 6.0}),
                space.makeVector(new double[]{5.0 / 6.0, 0, 3.0 / 6.0}),
                space.makeVector(new double[]{2.0 / 6.0, 0, 4.0 / 6.0}),
                space.makeVector(new double[]{5.0 / 6.0, 0, 5.0 / 6.0})});
        config.setBasisMolecules(alkaneBasis);
        config.setBasisSurface(surfaceBasis);
        config.setCellSizeX(sizeCellX);
        config.setCellSizeZ(sizeCellZ);
        config.setNCellsX(numXCells);
        config.setNCellsZ(numZCells);
        config.setSurfaceYOffset(2);

        config.setConformation(0, (ConformationChainZigZag) species.getConformation());
        config.setConformation(1, new ConformationChainZigZag(space));
        config.setConformation(2, new ConformationChainZigZag(space));
        config.setConformation(3, new ConformationChainZigZag(space));
        updateConformation(0);
        updateConformation(1);
        updateConformation(2);
        updateConformation(3);

        config.initializeCoordinates(box);

        integrator = new IntegratorVelocityVerletSAM(potentialMaster, random, 0.002, Kelvin.UNIT.toSim(300), box);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(500);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        AtomType typeCH2 = species.getCH2Type();
        AtomType typeCH3 = species.getCH3Type();
        AtomType typeS = species.getSulfurType();

        double sigmaCH3 = 3.75;
        double sigmaSulfur = 3.62;


        double epsilonCH2 = Kelvin.UNIT.toSim(46);
        double epsilonCH3 = Kelvin.UNIT.toSim(98);
        double epsilonSulfur = Kelvin.UNIT.toSim(232);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2 * epsilonCH3);
        double epsilonCH2Sulfur = Math.sqrt(epsilonCH2 * epsilonSulfur);
        // sulfur and CH3 will never be close
        double rc = 2.5 * sigmaCH2;
        double nbrCut = 2.8 * sigmaCH2;
        if (0.495 * box.getBoundary().getBoxSize().getX(0) < rc) {
            rc = 0.495 * box.getBoundary().getBoxSize().getX(0);
            nbrCut = 0.5 * box.getBoundary().getBoxSize().getX(0);
        }
        if (0.495 * box.getBoundary().getBoxSize().getX(2) < rc) {
            rc = 0.495 * box.getBoundary().getBoxSize().getX(2);
            nbrCut = 0.5 * box.getBoundary().getBoxSize().getX(2);
        }
        final double rCut = rc;
        potentialMaster.setRange(nbrCut);
        p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        p2S = new P2LennardJones(space, sigmaSulfur, epsilonSulfur);
        p2CH2CH3 = new P2LennardJones(space, 0.5 * (sigmaCH2 + sigmaCH3), epsilonCH2CH3);
        p2SCH2 = new P2LennardJones(space, 0.5 * (sigmaSulfur + sigmaCH2), epsilonCH2Sulfur);
        p2CH2t = new P2SoftSphericalTruncatedSwitched(space, p2CH2, rCut);
        p2CH3t = new P2SoftSphericalTruncatedSwitched(space, p2CH2, rCut);
        p2CH2CH3t = new P2SoftSphericalTruncatedSwitched(space, p2CH2, rCut);
        p2St = new P2SoftSphericalTruncatedSwitched(space, p2S, rCut);
        p2SCH2t = new P2SoftSphericalTruncatedSwitched(space, p2SCH2, rCut);

        NeighborCriterion nonBondedCriterion = new NeighborCriterion() {
            public boolean accept(IAtomList pair) {
                int idx0 = pair.get(0).getIndex();
                int idx1 = pair.get(1).getIndex();
                int idxDiff = idx0 - idx1;
                return idxDiff > 3 || idxDiff < -3;
            }

            public boolean needUpdate(IAtom atom) {
                return false;
            }

            public void reset(IAtom atom) {
            }

            public void setBox(Box box) {
            }

            public boolean unsafe() {
                return false;
            }
        };
        potentialMaster.addPotential(p2CH2t, new AtomType[]{typeCH2, typeCH2});
        ((CriterionInterMolecular) potentialMaster.getCriterion(typeCH2, typeCH2)).setIntraMolecularCriterion(nonBondedCriterion);
        potentialMaster.addPotential(p2CH3t, new AtomType[]{typeCH3, typeCH3});
        potentialMaster.addPotential(p2St, new AtomType[]{typeS, typeS});
        potentialMaster.addPotential(p2SCH2t, new AtomType[]{typeS, typeCH2});
        ((CriterionInterMolecular) potentialMaster.getCriterion(typeS, typeCH2)).setIntraMolecularCriterion(nonBondedCriterion);
        potentialMaster.addPotential(p2CH2CH3t, new AtomType[]{typeCH2, typeCH3});
        ((CriterionInterMolecular) potentialMaster.getCriterion(typeCH2, typeCH3)).setIntraMolecularCriterion(nonBondedCriterion);
        p1Intra = potentialMaster.makePotentialGroup(1);
        potentialMaster.addPotential(p1Intra, new ISpecies[]{species});

        p2BondCC = new P2Harmonic(space, 10000, bondL_CC);
        p2BondCS = new P2Harmonic(space, 10000, bondL_CS);
        // bond angle potential is the same for CCC and CCS
        p3Bond = new P3BondAngle(space);
        p3Bond.setAngle(Math.PI * 114.0 / 180.0);
        p3Bond.setEpsilon(Kelvin.UNIT.toSim(62500));
        p4BondCCCC = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
        p4BondCCCS = new P4BondTorsion(space, Kelvin.UNIT.toSim(-251.06), Kelvin.UNIT.toSim(428.73), Kelvin.UNIT.toSim(-111.85), Kelvin.UNIT.toSim(441.27));
        setChainLength(chainLength);
//        updateConformation(0);

        sinusoidalB = Calorie.UNIT.toSim(2000) / Constants.AVOGADRO;
        p1SurfaceBond = new P1Sinusoidal(space);
        p1SurfaceBond.setB(0); // initially disabled
        p1SurfaceBond.setCellSize(sizeCellX, sizeCellZ);
        p1SurfaceBond.setOffset(space.makeVector(new double[]{sizeCellX / 6.0, 0, sizeCellZ / 6.0}));
        potentialMaster.addPotential(p1SurfaceBond, new AtomType[]{species.getSulfurType()});

        wallPotential = new P1WCAWall(space, 1, 4, 1000);
        wallPotential.setWallPosition(box.getBoundary().getBoxSize().getX(1) * 0.5);
        potentialMaster.addPotential(wallPotential, new AtomType[]{species.getCH2Type()});
        potentialMaster.addPotential(wallPotential, new AtomType[]{species.getCH3Type()});

        forceSum = new PotentialCalculationForceSumWall(wallPotential);
        integrator.setForceSum(forceSum);

        P2LennardJones p2Surface = new P2LennardJones(space, 3.0, Kelvin.UNIT.toSim(50));
        p2CH2Surface = new P2SoftSphericalTruncatedSwitched(space, p2Surface, rCut);
        potentialMaster.addPotential(p2CH2Surface, new AtomType[]{speciesSurface.getLeafType(), species.getCH2Type()});
        potentialMaster.getNeighborManager(box).setDoApplyPBC(false);
        potentialMaster.getNbrCellManager(box).setDoApplyPBC(true);

        harmonicStrength = 10000;
        p2SurfaceBond = new P2Harmonic(space, harmonicStrength, 2.5);
        p2SulfurSurfaceLJ = new P2SoftSphericalTruncatedSwitched(space, p2Surface, rCut);
        potentialMaster.addPotential(p2SulfurSurfaceLJ, new AtomType[]{speciesSurface.getLeafType(), species.getSulfurType()});
        criterion3 = new CriterionTether3(this, species, speciesSurface.getLeafType());
        criterion3.setBox(box);
        P2Surface p2SurfaceTrunc = new P2Surface(space, p2SulfurSurfaceLJ, p2SurfaceBond, criterion3);
        potentialMaster.addPotential(p2SurfaceTrunc, new AtomType[]{speciesSurface.getLeafType(), species.getSulfurType()});
        findTetherBonds();


        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        updateRCut();
    }

    public static void main(String[] args) {
        Sam sim = new Sam();
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));
//        sim.integrator.setActionInterval(simGraphic.getPaintAction(sim.box), 10);
        simGraphic.getDisplayBox(sim.box).canvas.setDrawBoundary(DisplayCanvas.DRAW_BOUNDARY_NONE);
        simGraphic.makeAndDisplayFrame();
    }
    
    protected void updateConformation(int iChain) {
        double bondTheta0 = chainTheta + .5*(Math.PI - bondTheta);
        Vector vector1 = config.getConformation(iChain).getFirstVector();
        vector1.setX(0, Math.cos(chainPsi)*Math.sin(bondTheta0)*bondL_CC);
        vector1.setX(1, Math.cos(bondTheta0)*bondL_CC);
        vector1.setX(2, Math.sin(chainPsi)*Math.sin(bondTheta0)*bondL_CC);
        double bondTheta2 = bondTheta0 - (Math.PI - bondTheta);
        Vector vector2 = config.getConformation(iChain).getSecondVector();
        vector2.setX(0, Math.cos(chainPsi)*(Math.sin(bondTheta2))*bondL_CC);
        vector2.setX(1, Math.cos(bondTheta2)*bondL_CC);
        vector2.setX(2, Math.sin(chainPsi)*(Math.sin(bondTheta2))*bondL_CC);

        Vector vector0 = space.makeVector();
        vector0.Ev1Pv2(vector1, vector2);
        IOrientation3D orientation = (IOrientation3D)space.makeOrientation();
        orientation.setDirection(vector1);
        Vector vector0Axis = space.makeVector();
        vector0Axis.Ea1Tv1(1.0/Math.sqrt(vector0.squared()), vector0);
        orientation.rotateBy(chainPhi[iChain], vector0Axis);
        vector1.Ea1Tv1(Math.sqrt(vector1.squared()), orientation.getDirection());
        vector2.Ev1Mv2(vector0, vector1);

        if (iChain == 0) {
            IMolecule molecule = species.makeMolecule();
            Vector moleculePos = space.makeVector();
            moleculePos.E(positionDefinition.position(molecule));
            Vector sulfurPosition = molecule.getChildList().get(0).getPosition();
            sulfurPosition.ME(moleculePos);
            molecule = null;
            sulfurPosition.TE(-1);
            sulfurPosition.setX(1, sulfurPosition.getX(1)+2.5);

            config.setMoleculeOffset(sulfurPosition);
        }
    }

    public double getChainTheta() {
        return chainTheta;
    }
    
    public void setChainTheta(double newTheta) {
        chainTheta = newTheta;
        updateConformation(0);
        updateConformation(1);
        updateConformation(2);
        updateConformation(3);
    }

    public double getChainPsi() {
        return chainPsi;
    }
    
    public void setChainPsi(double newPsi) {
        chainPsi = newPsi;
        updateConformation(0);
        updateConformation(1);
        updateConformation(2);
        updateConformation(3);
    }
    
    public void setChainPhi(int iChain, double newPhi) {
        chainPhi[iChain] = newPhi;
        updateConformation(iChain);
    }
    
    public double getChainPhi(int iChain) {
        return chainPhi[iChain];
    }
    
    public void findTetherBonds() {
        IMoleculeList polymerMolecules = box.getMoleculeList(species);
        IMoleculeList surfaceMolecules = box.getMoleculeList(speciesSurface);
        int nMolecules = polymerMolecules.size();
        double maxDistance = 3.5*3.5;
        Vector dr = space.makeVector();
        Boundary boundary = box.getBoundary();
        for (int i=0; i<nMolecules; i++) {
            AtomArrayList bondedSurfaceAtoms = new AtomArrayList(3);
            IAtom sulfur = polymerMolecules.get(i).getChildList().get(0);
            for (int j = 0; j<surfaceMolecules.size(); j++) {
                IAtom gold = surfaceMolecules.get(j).getChildList().get(0);
                dr.Ev1Mv2(sulfur.getPosition(), gold.getPosition());
                boundary.nearestImage(dr);
                if (dr.squared() < maxDistance) {
                    bondedSurfaceAtoms.add(gold);
                }
            }
            if (bondedSurfaceAtoms.size() != 3) {
                throw new RuntimeException("only found "+bondedSurfaceAtoms.size()+" bonded atoms");
            }
            criterion3.setBondedSurfaceAtoms(polymerMolecules.get(i), bondedSurfaceAtoms);
        }
    }

    public int getNumZCells() {
        return numZCells;
    }

    public void setNumZCells(int newNumZCells) {
        if (newNumZCells == numZCells) {
            return;
        }
        boolean increase = newNumZCells > numZCells;
        numZCells = newNumZCells;
        Vector dim = space.makeVector();
        double zShift = box.getBoundary().getBoxSize().getX(2);
        dim.E(new double[]{sizeCellX*numXCells, chainLength*2.5, sizeCellZ*numZCells});
        box.getBoundary().setBoxSize(dim);
        config.setNCellsZ(numZCells);
        if (!increase) {
            config.initializeCoordinates(box);
        }
        else {
            IAtomList leafList = box.getLeafList();
            for (int i = 0; i<leafList.size(); i++) {
                IAtom a = leafList.get(i);
                a.getPosition().setX(2, a.getPosition().getX(2) - 0.5*zShift);
            }

            box.setNMolecules(species, 2*numXCells*numZCells);
            box.setNMolecules(speciesSurface, 6*numXCells*numZCells);
            IMoleculeList molecules = box.getMoleculeList(species);
            for (int i = 0; i<molecules.size()/2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtomList childList = molecules.get(i+molecules.size()/2).getChildList();
                for (int j = 0; j<childList.size(); j++) {
                    IAtom atom0 = childList0.get(j);
                    IAtom atom = childList.get(j);
                    atom.getPosition().E(atom0.getPosition());
                    atom.getPosition().setX(2, atom.getPosition().getX(2) + zShift);
                }
            }

            molecules = box.getMoleculeList(speciesSurface);
            for (int i = 0; i<molecules.size()/2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtom atom0 = childList0.get(0);
                IAtomList childList = molecules.get(i+molecules.size()/2).getChildList();
                IAtom atom = childList.get(0);
                atom.getPosition().E(atom0.getPosition());
                atom.getPosition().setX(2, atom.getPosition().getX(2) + zShift);
            }
        }

        updateRCut();
        findTetherBonds();

        integrator.reset();
    }

    public int getNumXCells() {
        return numXCells;
    }

    public void setNumXCells(int newNumXCells) {
        if (newNumXCells == numXCells) {
            return;
        }
        int oldNumXCells = numXCells;
        numXCells = newNumXCells;
        double xShift = box.getBoundary().getBoxSize().getX(0);
        Vector dim = space.makeVector();
        dim.E(new double[]{sizeCellX*numXCells, chainLength*2.5, sizeCellZ*numZCells});
        box.getBoundary().setBoxSize(dim);
        config.setNCellsX(numXCells);
        if (numXCells == 2) {
            p1SurfaceBond.setOffset(space.makeVector(new double[]{4*sizeCellX/6.0, 0, sizeCellZ/6.0}));
        }
        else {
            p1SurfaceBond.setOffset(space.makeVector(new double[]{sizeCellX/6.0, 0, sizeCellZ/6.0}));
        }
        if (numXCells < oldNumXCells || oldNumXCells*2 != numXCells) {
            config.initializeCoordinates(box);
        }
        else {
            IAtomList leafList = box.getLeafList();
            for (int i = 0; i<leafList.size(); i++) {
                IAtom a = leafList.get(i);
                a.getPosition().setX(0, a.getPosition().getX(0) - 0.5*xShift);
            }

            box.setNMolecules(species, 2*numXCells*numZCells);
            box.setNMolecules(speciesSurface, 6*numXCells*numZCells);
            IMoleculeList molecules = box.getMoleculeList(species);
            for (int i = 0; i<molecules.size()/2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtomList childList = molecules.get(i+molecules.size()/2).getChildList();
                for (int j = 0; j<childList.size(); j++) {
                    IAtom atom0 = childList0.get(j);
                    IAtom atom = childList.get(j);
                    atom.getPosition().E(atom0.getPosition());
                    atom.getPosition().setX(0, atom.getPosition().getX(0) + xShift);
                }
            }

            molecules = box.getMoleculeList(speciesSurface);
            for (int i = 0; i<molecules.size()/2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtom atom0 = childList0.get(0);
                IAtomList childList = molecules.get(i+molecules.size()/2).getChildList();
                IAtom atom = childList.get(0);
                atom.getPosition().E(atom0.getPosition());
                atom.getPosition().setX(0, atom.getPosition().getX(0) + xShift);
            }
        }

        updateRCut();
        findTetherBonds();
        integrator.reset();
    }

    protected void updateRCut() {
        double rCut = 2.5*sigmaCH2;
        double nbrCut = 2.8*sigmaCH2;
        if (0.5*box.getBoundary().getBoxSize().getX(0) < nbrCut) {
            nbrCut = 0.5*box.getBoundary().getBoxSize().getX(0);
            rCut = nbrCut * 0.85;
        }
        if (0.5*box.getBoundary().getBoxSize().getX(2) < nbrCut) {
            nbrCut = 0.5*box.getBoundary().getBoxSize().getX(2);
            rCut = nbrCut * 0.85;
        }
        p2CH2t.setTruncationRadius(rCut);
        p2CH3t.setTruncationRadius(rCut);
        p2CH2CH3t.setTruncationRadius(rCut);
        p2St.setTruncationRadius(rCut);
        p2SCH2t.setTruncationRadius(rCut);
        p2SulfurSurfaceLJ.setTruncationRadius(rCut);
        p2CH2Surface.setTruncationRadius(rCut);

        potentialMaster.setRange(nbrCut);
        potentialMaster.reset();
    }

    public void setChainLength(int newChainLength) {
        if (newChainLength < 7) {
            throw new RuntimeException("too short!");
        }
        chainLength = newChainLength;
        potentialMaster.removePotential(p1Intra);
        p1Intra = integrator.getPotentialMaster().makePotentialGroup(1);
        potentialMaster.addPotential(p1Intra, new ISpecies[]{species});

        int[][] pairs = new int[1][2];
        pairs[0][0] = 0;
        pairs[0][1] = 1;
        p1Intra.addPotential(p2BondCS, new ApiIndexList(pairs));
        pairs = new int[chainLength-2][2];
        for (int i=1; i<chainLength-1; i++) {
            pairs[i-1][0] = i;
            pairs[i-1][1] = i+1;
        }
        p1Intra.addPotential(p2BondCC, new ApiIndexList(pairs));

        // CCC and CCS bond is the same
        int[][] triplets = new int[chainLength-2][3];
        for (int i=0; i<chainLength-2; i++) {
            triplets[i][0] = i;
            triplets[i][1] = i+1;
            triplets[i][2] = i+2;
        }
        p1Intra.addPotential(p3Bond, new Atomset3IteratorIndexList(triplets));

        int[][] quads = new int[][]{{0,1,2,3}};
        p1Intra.addPotential(p4BondCCCS, new Atomset4IteratorIndexList(quads));
        quads = new int[chainLength-3][4];
        for (int i=0; i<chainLength-3; i++) {
            quads[i][0] = i;
            quads[i][1] = i+1;
            quads[i][2] = i+2;
            quads[i][3] = i+3;
        }
        p1Intra.addPotential(p4BondCCCC, new Atomset4IteratorIndexList(quads));
    }

    /**
     * Potential class that behaves like switched LJ except for the specific pair
     * of sulfur/surface atoms that are tethered.
     */
    private static class P2Surface extends Potential2 implements PotentialSoft {
        protected final PotentialSoft p2lj;
        protected final P2Harmonic p2Bond;
        protected final Vector dr;
        protected Boundary boundary;
        protected final NeighborCriterion bondCriterion;

        public P2Surface(Space space, PotentialSoft p2lj, P2Harmonic p2Bond, NeighborCriterion bondCriterion) {
            super(space);
            dr = space.makeVector();
            this.p2lj = p2lj;
            this.p2Bond = p2Bond;
            this.bondCriterion = bondCriterion;
        }

        public double getRange() {
            return p2lj.getRange();
        }

        @Override
        public double energy(IAtomList atoms) {
            return bondCriterion.accept(atoms) ? p2Bond.energy(atoms) : p2lj.energy(atoms);
        }

        @Override
        public void setBox(Box box) {
            p2lj.setBox(box);
            p2Bond.setBox(box);
        }

        @Override
        public double virial(IAtomList atoms) {
            return bondCriterion.accept(atoms) ? p2Bond.virial(atoms) : p2lj.virial(atoms);
        }

        @Override
        public Vector[] gradient(IAtomList atoms) {
            return bondCriterion.accept(atoms) ? p2Bond.gradient(atoms) : p2lj.gradient(atoms);
        }

        @Override
        public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
            return bondCriterion.accept(atoms) ? p2Bond.gradient(atoms, pressureTensor) : p2lj.gradient(atoms, pressureTensor);
        }
    }
}
