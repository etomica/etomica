/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Sulfur;
import etomica.config.ConformationChainZigZag;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.lattice.crystal.Basis;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManagerFasterer;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.*;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.IOrientation3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Calorie;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.Constants;

import java.util.ArrayList;
import java.util.List;

/**
 * Self-assembled monolayer module.
 *
 * @author Andrew Schultz
 */
public class SamFasterer extends Simulation {

    public final P1Sinusoidal p1SurfaceBond;
    public final double sinusoidalB;
    public final P2Harmonic p2SurfaceBond;
    public final double harmonicStrength;
    public SpeciesGeneral species;
    public SpeciesGeneral speciesSurface;
    public Box box;
    public IntegratorVelocityVerletFasterer integrator;

    public IMoleculePositionDefinition positionDefinition;
    public P1WCAWall wallPotential;
    public ConfigurationSAM config;
    public P2LennardJones p2CH2, p2CH3, p2CH2CH3, p2S, p2SCH2;
    public P2SoftSphericalTruncatedSwitched p2CH2t, p2CH3t, p2CH2CH3t, p2St, p2SCH2t;
    public P2Harmonic p2BondCC, p2BondCS;
    public P3BondAngle p3Bond;
    public P4BondTorsion p4BondCCCC, p4BondCCCS;
    public P2SoftSphericalTruncatedSwitched p2SulfurSurfaceLJ, p2CH2Surface;
    public CriterionTether3 criterion3;
    public int chainLength;
    public double chainTheta, chainPsi;
    public double[] chainPhi;
    public double bondL_CC;
    public double bondTheta;
    public PotentialCalculationForceSumWall forceSum;
    public PotentialMasterListFasterer potentialMaster;
    public PotentialMasterBonding potentialMasterBonding;
    public PotentialComputePairGeneral potentialMasterGeneral;
    public NeighborListManagerFasterer neighborManagerGeneral;
    public int numZCells, numXCells;
    public double sizeCellZ, sizeCellX;
    public double sigmaCH2;
    public SamIntegratorListener integratorListener;
    public double wallForce;

    public SamFasterer() {
        super(Space.getInstance(3));

        //species and potentials
        chainLength = 16;
        ConformationChainZigZag conf0 = new ConformationChainZigZag(getSpace());
        species = new SpeciesBuilder(getSpace())
                .addCount(new AtomType(Sulfur.INSTANCE), 1)
                .addCount(AtomType.simple("CH2", 14), chainLength - 2)
                .addCount(AtomType.simple("CH3", 15), 1)
                .setDynamic(true)
                .withConformation(conf0)
                .build();
        addSpecies(species);

        speciesSurface = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) speciesSurface.getLeafType().getElement()).setMass(Double.POSITIVE_INFINITY);
        addSpecies(speciesSurface);

        //construct box
        box = this.makeBox(new BoundaryRectangularSlit(1, space));
        Vector dim = space.makeVector();
        dim.E(new double[]{sizeCellX * numXCells, chainLength * 2.4, sizeCellZ * numZCells});
        box.getBoundary().setBoxSize(dim);

        PotentialComputeAggregate pcAggregate = new PotentialComputeAggregate();
        potentialMasterBonding = new PotentialMasterBonding(this, box);
        pcAggregate.add(potentialMasterBonding);

        positionDefinition = new MoleculePositionGeometricCenter(space);
        sigmaCH2 = 3.95;
        potentialMaster = new PotentialMasterListFasterer(getSpeciesManager(), box, 2, 2.8 * sigmaCH2, potentialMasterBonding.getBondingInfo()); //List(this, 2.0);
        pcAggregate.add(potentialMaster);

        numXCells = 4;
        numZCells = 2;
        // gold has FCC unit cell, a=4.0782A
        sizeCellZ = 4.0782 / Math.sqrt(2) * 3; //Math.sqrt(3)*sizeCellX;
        sizeCellX = sizeCellZ / Math.sqrt(3);

        bondL_CC = 1.54;
        double bondL_CS = 1.82;
        bondTheta = Math.PI * 114 / 180;
        chainTheta = 0;
        chainPsi = 0;
        chainPhi = new double[4];

        config = new ConfigurationSAM(this, space, species, speciesSurface, null);
        Basis alkaneBasis = new Basis(new Vector[]{Vector.of(new double[]{1.0 / 6.0, 0, 1.0 / 6.0}), Vector.of(new double[]{2.0 / 3.0, 0, 2.0 / 3.0})});
        Basis surfaceBasis = new Basis(new Vector[]{
                Vector.of(new double[]{2.0 / 6.0, 0, 0}),
                Vector.of(new double[]{5.0 / 6.0, 0, 1.0 / 6.0}),
                Vector.of(new double[]{2.0 / 6.0, 0, 2.0 / 6.0}),
                Vector.of(new double[]{5.0 / 6.0, 0, 3.0 / 6.0}),
                Vector.of(new double[]{2.0 / 6.0, 0, 4.0 / 6.0}),
                Vector.of(new double[]{5.0 / 6.0, 0, 5.0 / 6.0})});
        config.setBasisMolecules(alkaneBasis);
        config.setBasisSurface(surfaceBasis);
        config.setCellSizeX(sizeCellX);
        config.setCellSizeZ(sizeCellZ);
        config.setNCellsX(numXCells);
        config.setNCellsZ(numZCells);
        config.setSurfaceYOffset(2);

        config.setConformation(0, conf0);
        config.setConformation(1, new ConformationChainZigZag(space));
        config.setConformation(2, new ConformationChainZigZag(space));
        config.setConformation(3, new ConformationChainZigZag(space));
        updateConformation(0);
        updateConformation(1);
        updateConformation(2);
        updateConformation(3);

        config.initializeCoordinates(box);

        AtomType typeCH2 = species.getTypeByName("CH2");
        AtomType typeCH3 = species.getTypeByName("CH3");
        AtomType typeS = species.getTypeByName("S");

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
        potentialMaster.setNeighborRange(nbrCut);
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

        potentialMaster.setPairPotential(typeCH2, typeCH2, p2CH2t);
        potentialMaster.setPairPotential(typeCH3, typeCH3, p2CH3t);
        potentialMaster.setPairPotential(typeS, typeS, p2St);
        potentialMaster.setPairPotential(typeS, typeCH2, p2SCH2t);
        potentialMaster.setPairPotential(typeCH2, typeCH3, p2CH2CH3t);

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
        p1SurfaceBond.setOffset(Vector.of(new double[]{sizeCellX / 6.0, 0, sizeCellZ / 6.0}));
        PotentialComputeField computeField = new PotentialComputeField(this, box) {
            public double computeAll(boolean doForces, PotentialCallback pc) {
                double u = super.computeAll(doForces, pc);
                wallForce = 0;
                for (Vector f : forces) {
                    wallForce += f.getX(1);
                }
                return u;
            }
        };
        computeField.setFieldPotential(typeS, p1SurfaceBond);

        wallPotential = new P1WCAWall(space, 1, 4, 1000);
        wallPotential.setWallPosition(box.getBoundary().getBoxSize().getX(1) * 0.5);
        computeField.setFieldPotential(typeCH2, wallPotential);
        computeField.setFieldPotential(typeCH3, wallPotential);
        pcAggregate.add(computeField);

        P2LennardJones p2Surface = new P2LennardJones(space, 3.0, Kelvin.UNIT.toSim(50));
        p2CH2Surface = new P2SoftSphericalTruncatedSwitched(space, p2Surface, rCut);
        potentialMaster.setPairPotential(speciesSurface.getLeafType(), typeCH2, p2CH2Surface);

        harmonicStrength = 10000;
        p2SurfaceBond = new P2Harmonic(space, harmonicStrength, 2.5);
        p2SulfurSurfaceLJ = new P2SoftSphericalTruncatedSwitched(space, p2Surface, rCut);
        potentialMaster.setPairPotential(speciesSurface.getLeafType(), typeS, p2SulfurSurfaceLJ);
        criterion3 = new CriterionTether3(this, species, speciesSurface.getLeafType());
        criterion3.setBox(box);
        P2Surface p2SurfaceTrunc = new P2Surface(space, p2SulfurSurfaceLJ, p2SurfaceBond, criterion3);

        neighborManagerGeneral = new NeighborListManagerFasterer(getSpeciesManager(), box, 2, nbrCut, potentialMasterBonding.getBondingInfo());
        potentialMasterGeneral = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManagerGeneral);
        potentialMasterGeneral.setPairPotential(speciesSurface.getLeafType(), typeS, p2SurfaceTrunc);
        pcAggregate.add(potentialMasterGeneral);
        findTetherBonds();

        integrator = new IntegratorVelocityVerletFasterer(pcAggregate, random, 0.002, Kelvin.UNIT.toSim(300), box);

        integrator.setIsothermal(true);
        integrator.setThermostatInterval(500);
        getController().addActivity(new ActivityIntegrate(integrator));

        forceSum = new PotentialCalculationForceSumWall(wallPotential);
        integratorListener = new SamIntegratorListener(pcAggregate, box);
        integrator.getEventManager().addListener(integratorListener);

        updateRCut();
    }

    public static void main(String[] args) {
        SamFasterer sim = new SamFasterer();
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));
//        sim.integrator.setActionInterval(simGraphic.getPaintAction(sim.box), 10);
        simGraphic.getDisplayBox(sim.box).canvas.setDrawBoundary(DisplayCanvas.DRAW_BOUNDARY_NONE);
        simGraphic.makeAndDisplayFrame();
    }

    protected void updateConformation(int iChain) {
        double bondTheta0 = chainTheta + .5 * (Math.PI - bondTheta);
        Vector vector1 = config.getConformation(iChain).getFirstVector();
        vector1.setX(0, Math.cos(chainPsi) * Math.sin(bondTheta0) * bondL_CC);
        vector1.setX(1, Math.cos(bondTheta0) * bondL_CC);
        vector1.setX(2, Math.sin(chainPsi) * Math.sin(bondTheta0) * bondL_CC);
        double bondTheta2 = bondTheta0 - (Math.PI - bondTheta);
        Vector vector2 = config.getConformation(iChain).getSecondVector();
        vector2.setX(0, Math.cos(chainPsi) * (Math.sin(bondTheta2)) * bondL_CC);
        vector2.setX(1, Math.cos(bondTheta2) * bondL_CC);
        vector2.setX(2, Math.sin(chainPsi) * (Math.sin(bondTheta2)) * bondL_CC);

        Vector vector0 = space.makeVector();
        vector0.Ev1Pv2(vector1, vector2);
        IOrientation3D orientation = (IOrientation3D) space.makeOrientation();
        orientation.setDirection(vector1);
        Vector vector0Axis = space.makeVector();
        vector0Axis.Ea1Tv1(1.0 / Math.sqrt(vector0.squared()), vector0);
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
            sulfurPosition.setX(1, sulfurPosition.getX(1) + 2.5);

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
        double maxDistance = 3.5 * 3.5;
        Vector dr = space.makeVector();
        Boundary boundary = box.getBoundary();
        for (int i = 0; i < nMolecules; i++) {
            AtomArrayList bondedSurfaceAtoms = new AtomArrayList(3);
            IAtom sulfur = polymerMolecules.get(i).getChildList().get(0);
            for (int j = 0; j < surfaceMolecules.size(); j++) {
                IAtom gold = surfaceMolecules.get(j).getChildList().get(0);
                dr.Ev1Mv2(sulfur.getPosition(), gold.getPosition());
                boundary.nearestImage(dr);
                if (dr.squared() < maxDistance) {
                    bondedSurfaceAtoms.add(gold);
                }
            }
            if (bondedSurfaceAtoms.size() != 3) {
                throw new RuntimeException("only found " + bondedSurfaceAtoms.size() + " bonded atoms");
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
        dim.E(new double[]{sizeCellX * numXCells, chainLength * 2.5, sizeCellZ * numZCells});
        box.getBoundary().setBoxSize(dim);
        config.setNCellsZ(numZCells);
        if (!increase) {
            config.initializeCoordinates(box);
        } else {
            IAtomList leafList = box.getLeafList();
            for (int i = 0; i < leafList.size(); i++) {
                IAtom a = leafList.get(i);
                a.getPosition().setX(2, a.getPosition().getX(2) - 0.5 * zShift);
            }

            box.setNMolecules(species, 2 * numXCells * numZCells);
            box.setNMolecules(speciesSurface, 6 * numXCells * numZCells);
            IMoleculeList molecules = box.getMoleculeList(species);
            for (int i = 0; i < molecules.size() / 2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtomList childList = molecules.get(i + molecules.size() / 2).getChildList();
                for (int j = 0; j < childList.size(); j++) {
                    IAtom atom0 = childList0.get(j);
                    IAtom atom = childList.get(j);
                    atom.getPosition().E(atom0.getPosition());
                    atom.getPosition().setX(2, atom.getPosition().getX(2) + zShift);
                }
            }

            molecules = box.getMoleculeList(speciesSurface);
            for (int i = 0; i < molecules.size() / 2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtom atom0 = childList0.get(0);
                IAtomList childList = molecules.get(i + molecules.size() / 2).getChildList();
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
        dim.E(new double[]{sizeCellX * numXCells, chainLength * 2.5, sizeCellZ * numZCells});
        box.getBoundary().setBoxSize(dim);
        config.setNCellsX(numXCells);
        if (numXCells == 2) {
            p1SurfaceBond.setOffset(Vector.of(new double[]{4 * sizeCellX / 6.0, 0, sizeCellZ / 6.0}));
        } else {
            p1SurfaceBond.setOffset(Vector.of(new double[]{sizeCellX / 6.0, 0, sizeCellZ / 6.0}));
        }
        if (numXCells < oldNumXCells || oldNumXCells * 2 != numXCells) {
            config.initializeCoordinates(box);
        } else {
            IAtomList leafList = box.getLeafList();
            for (int i = 0; i < leafList.size(); i++) {
                IAtom a = leafList.get(i);
                a.getPosition().setX(0, a.getPosition().getX(0) - 0.5 * xShift);
            }

            box.setNMolecules(species, 2 * numXCells * numZCells);
            box.setNMolecules(speciesSurface, 6 * numXCells * numZCells);
            IMoleculeList molecules = box.getMoleculeList(species);
            for (int i = 0; i < molecules.size() / 2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtomList childList = molecules.get(i + molecules.size() / 2).getChildList();
                for (int j = 0; j < childList.size(); j++) {
                    IAtom atom0 = childList0.get(j);
                    IAtom atom = childList.get(j);
                    atom.getPosition().E(atom0.getPosition());
                    atom.getPosition().setX(0, atom.getPosition().getX(0) + xShift);
                }
            }

            molecules = box.getMoleculeList(speciesSurface);
            for (int i = 0; i < molecules.size() / 2; i++) {
                IAtomList childList0 = molecules.get(i).getChildList();
                IAtom atom0 = childList0.get(0);
                IAtomList childList = molecules.get(i + molecules.size() / 2).getChildList();
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
        double rCut = 2.5 * sigmaCH2;
        double nbrCut = 2.8 * sigmaCH2;
        if (0.5 * box.getBoundary().getBoxSize().getX(0) < nbrCut) {
            nbrCut = 0.5 * box.getBoundary().getBoxSize().getX(0);
            rCut = nbrCut * 0.85;
        }
        if (0.5 * box.getBoundary().getBoxSize().getX(2) < nbrCut) {
            nbrCut = 0.5 * box.getBoundary().getBoxSize().getX(2);
            rCut = nbrCut * 0.85;
        }
        p2CH2t.setTruncationRadius(rCut);
        p2CH3t.setTruncationRadius(rCut);
        p2CH2CH3t.setTruncationRadius(rCut);
        p2St.setTruncationRadius(rCut);
        p2SCH2t.setTruncationRadius(rCut);
        p2SulfurSurfaceLJ.setTruncationRadius(rCut);
        p2CH2Surface.setTruncationRadius(rCut);

        potentialMaster.setNeighborRange(nbrCut);
        neighborManagerGeneral.setNeighborRange(nbrCut);
        potentialMaster.init();
    }

    public void setChainLength(int newChainLength) {
        if (newChainLength < 7) {
            throw new RuntimeException("too short!");
        }
        chainLength = newChainLength;

        List<int[]> pairs = new ArrayList<>();
        pairs.add(new int[]{0, 1});
        potentialMasterBonding.setBondingPotentialPair(species, p2BondCS, pairs);
        pairs = new ArrayList<>();
        for (int i = 1; i < chainLength - 1; i++) {
            pairs.add(new int[]{i, i + 1});
        }
        potentialMasterBonding.setBondingPotentialPair(species, p2BondCC, pairs);

        // CCC and CCS bond is the same
        List<int[]> triplets = new ArrayList<>();
        for (int i = 0; i < chainLength - 2; i++) {
            triplets.add(new int[]{i, i + 1, i + 2});
        }
        potentialMasterBonding.setBondingPotentialTriplet(species, p3Bond, triplets);

        List<int[]> quads = new ArrayList<>();
        quads.add(new int[]{0, 1, 2, 3});
        potentialMasterBonding.setBondingPotentialQuad(species, p4BondCCCS, quads);
        quads = new ArrayList<>();
        for (int i = 0; i < chainLength - 3; i++) {
            quads.add(new int[]{i, i + 1, i + 2, i + 3});
        }
        potentialMasterBonding.setBondingPotentialQuad(species, p4BondCCCC, quads);
    }

    /**
     * Potential class that behaves like switched LJ except for the specific pair
     * of sulfur/surface atoms that are tethered.
     */
    private static class P2Surface extends Potential2 implements PotentialSoft, IPotentialPair {
        protected final Potential2Soft p2lj;
        protected final P2Harmonic p2Bond;
        protected final Vector dr;
        protected Boundary boundary;
        protected final NeighborCriterion bondCriterion;

        public P2Surface(Space space, Potential2Soft p2lj, P2Harmonic p2Bond, NeighborCriterion bondCriterion) {
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
            return bondCriterion.accept(atoms.get(0), atoms.get(1)) ? p2Bond.energy(atoms) : p2lj.energy(atoms);
        }

        @Override
        public void setBox(Box box) {
            p2lj.setBox(box);
            p2Bond.setBox(box);
        }

        @Override
        public double virial(IAtomList atoms) {
            return bondCriterion.accept(atoms.get(0), atoms.get(1)) ? p2Bond.virial(atoms) : p2lj.virial(atoms);
        }

        @Override
        public Vector[] gradient(IAtomList atoms) {
            return bondCriterion.accept(atoms.get(0), atoms.get(1)) ? p2Bond.gradient(atoms) : p2lj.gradient(atoms);
        }

        @Override
        public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
            return bondCriterion.accept(atoms.get(0), atoms.get(1)) ? p2Bond.gradient(atoms, pressureTensor) : p2lj.gradient(atoms, pressureTensor);
        }

        @Override
        public double u(Vector dr12, IAtom atom1, IAtom atom2) {
            double r2 = dr12.squared();
            return bondCriterion.accept(atom1, atom2) ? p2Bond.u(r2) : p2lj.u(r2);
        }

        @Override
        public double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
            double r2 = dr12.squared();
            double[] u012 = new double[3];
            if (bondCriterion.accept(atom1, atom2)) {
                p2Bond.u012add(r2, u012);
            } else {
                p2lj.u012add(r2, u012);
            }
            if (u012[1] == 0) return u012[0];
            f2.PEa1Tv1(-u012[1] / r2, dr12);
            f1.PEa1Tv1(u012[1] / r2, dr12);
            return u012[0];
        }

        @Override
        public Vector[][] gradientAndTorque(IAtomList atoms) {
            throw new RuntimeException("nope");
        }
    }
}
