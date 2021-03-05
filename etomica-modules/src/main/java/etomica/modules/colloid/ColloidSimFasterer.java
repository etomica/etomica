/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
import etomica.potential.IPotentialAtomic;
import etomica.potential.P1HardFieldGeneric;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;

/**
 * Colloid simulation.  Design by Alberto Striolo.
 *
 * @author Andrew Schultz
 */
public class ColloidSimFasterer extends Simulation {

    public PotentialComputePair potentialMaster;
    public SpeciesGeneral species, speciesColloid;
    public Box box;
    public IntegratorHardFasterer integrator;

    public ConfigurationColloid configuration;
    public AtomLeafAgentManager<AtomArrayList> colloidMonomerBondManager;
    public AtomLeafAgentManager<AtomArrayList> monomerMonomerBondManager;
    public P2SquareWellMonomerFasterer p2mm;
    public P2HardSphereMCFasterer p2mc;
    public int nGraft;
    public int chainLength;
    public P1WallFasterer p1WallMonomer;
    public P1HardFieldGeneric p1WallColloid;
    public double epsWallWall = 0.25;
    public double rangeWallColloid;

    public ColloidSimFasterer() {
        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(2));
        //species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        speciesColloid = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesColloid);

        box = this.makeBox();

        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, 6, BondingInfo.noBonding()) {
            protected int checkNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, double rc2, Vector jbo, IPotentialAtomic[] iPotentials) {
                int idx1 = iAtom.getParentGroup().getIndex();
                int idx2 = jAtom.getParentGroup().getIndex();
                int chainIdx1 = idx1 / chainLength;
                int chainIdx2 = idx2 / chainLength;
                int childIdx1 = idx1 % chainLength;
                int childIdx2 = idx2 % chainLength;
                boolean grafted = childIdx1 + childIdx2 == 0;
                boolean bonded = (chainIdx1 == chainIdx2) && Math.abs(childIdx1 - childIdx2) == 1;
                if (grafted || bonded) {
                    Vector dr = space.makeVector();
                    Vector ri = iAtom.getPosition();
                    Vector rj = jAtom.getPosition();
                    dr.Ev1Mv2(rj, ri);
                    dr.PE(jbo);
                    return addAsNbrPair(i, j, iAtom, jAtom, jbo, iPotentials, dr);
                }
                return super.checkNbrPair(i, j, iAtom, jAtom, rc2, jbo, iPotentials);
            }
        };
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        int nColloid = 1;
        chainLength = 50;
        nGraft = 12;

        double sigma = 1.0;
        double sigmaColloid = 7.5;
        double lambda = 1.5;
        double boxSize = 80;
        double epsMM = 1.0;

        p2mm = new P2SquareWellMonomerFasterer(sigma, lambda, epsMM, 0.85, chainLength);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2mm);
        p2mc = new P2HardSphereMCFasterer(0.5 * (sigma + sigmaColloid), 0.9, chainLength, colloidMonomerBondManager, species);
        potentialMaster.setPairPotential(species.getLeafType(), speciesColloid.getLeafType(), p2mc);

        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        double L = boxSize;
        double epsMWall = Math.sqrt(epsMM * epsWallWall);

        p1WallMonomer = new P1WallFasterer(box, new double[]{-L / 2 + 1.0 / 2, -L / 2 + 2,
                +L / 2 - 2, +L / 2 - 1.0 / 2}, new double[]{-epsMWall, 0, -epsMWall}, chainLength);
        pcField.setFieldPotential(species.getLeafType(), p1WallMonomer);

        rangeWallColloid = 10;
        p1WallColloid = new P1HardFieldGeneric(1, new double[]{-L / 2 + 7.5 / 2, -L / 2 + rangeWallColloid,
                +L / 2 - rangeWallColloid, +L / 2 - 7.5 / 2}, new double[]{-0.5, 0, -0.5}, false);
        pcField.setFieldPotential(speciesColloid.getLeafType(), p1WallColloid);


        //controller and integrator
        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), IntegratorHardFasterer.extractFieldPotentials(pcField), neighborManager, random, 0.02, 2, box, null, null);
        integrator.setTimeStep(0.02);
        integrator.setTemperature(2);
        integrator.setIsothermal(true);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);

        //potentials
        Vector dim = space.makeVector();
        dim.E(boxSize);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(speciesColloid, nColloid);

        AgentSource<AtomArrayList> bondAgentSource = new AgentSource<AtomArrayList>() {
            public void releaseAgent(AtomArrayList agent, IAtom atom, Box agentBox) {
            }

            public AtomArrayList makeAgent(IAtom a, Box agentBox) {
                return new AtomArrayList();
            }
        };
        colloidMonomerBondManager = new AtomLeafAgentManager<>(bondAgentSource, box);
        monomerMonomerBondManager = new AtomLeafAgentManager<>(bondAgentSource, box);

        //instantiate several potentials for selection in combo-box
//        potentialMaster.setCriterion(species.getLeafType(), speciesColloid.getLeafType(), new CriterionAll() {
//            public boolean needUpdate(IAtom atom) {
//                // just do this to force PBC to be applied
//                return integrator.getStepCount() % 100 == 0;
//            }
//        });


        //construct box
        configuration = new ConfigurationColloid(space, species, speciesColloid, random);
        configuration.setNGraft(nGraft);
        configuration.setChainLength(chainLength);
        configuration.setSigmaColloid(sigmaColloid);
        configuration.setSigmaMonomer(sigma);
        configuration.setMonomerMonomerBondManager(monomerMonomerBondManager);
        configuration.setColloidMonomerBondManager(colloidMonomerBondManager);
        configuration.initializeCoordinates(box);
        configuration.initializeCoordinates(box);

        if (nGraft > 1) {
            // the whole point here is to keep the chains reasonably distributed on the surface of the colloid
            // we add here an HS potential between the monomers that are grafted
            IAtom atom1 = box.getMoleculeList(species).get(0).getChildList().get(0);
            double minr2 = Double.POSITIVE_INFINITY;
            for (int j = chainLength; j < box.getMoleculeList(species).size(); j += chainLength) {
                IAtom atom2 = box.getMoleculeList(species).get(j).getChildList().get(0);
                Vector dr = space.makeVector();
                dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                if (minr2 > r2) {
                    minr2 = r2;
                }
            }
            double rMin = 0.9 * Math.sqrt(minr2);
            p2mm.setRGraftMin(rMin);
        }

    }

    public IntegratorBoxFasterer getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        ColloidSimFasterer sim = new ColloidSimFasterer();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }

    public void setNumGraft(int newNumGraft) {
        if (nGraft == newNumGraft) return;
        configuration.setNGraft(newNumGraft);
        nGraft = newNumGraft;
        configuration.initializeCoordinates(box);

        if (newNumGraft > 1) {

            IAtom atom1 = box.getMoleculeList(species).get(0).getChildList().get(0);
            double minr2 = Double.POSITIVE_INFINITY;
            for (int j = chainLength; j < box.getMoleculeList(species).size(); j += chainLength) {
                IAtom atom2 = box.getMoleculeList(species).get(j).getChildList().get(0);
                Vector dr = space.makeVector();
                dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                if (minr2 > r2) {
                    minr2 = r2;
                }
            }
            p2mm.setRGraftMin(0.9 * Math.sqrt(minr2));
        }
        try {
            integrator.reset();
        } catch (ConfigurationOverlapException e) {
        }
    }

    public int getChainLength() {
        return chainLength;
    }

    public void setChainLength(int newChainLength) {
        if (chainLength == newChainLength) return;
        configuration.setChainLength(newChainLength);
        chainLength = newChainLength;
        configuration.initializeCoordinates(box);
        p2mm.setChainLength(newChainLength);
        p2mc.setChainLength(newChainLength);
        p1WallMonomer.setChainLength(newChainLength);
//        ((NeighborListManagerColloid)potentialMaster.getNeighborManager(box)).setChainLength(chainLength);

        try {
            integrator.reset();
        } catch (ConfigurationOverlapException e) {
        }
    }

    public double getColloidSigma() {
        double coreMC = p2mc.getCollisionDiameter(0);
        double coreMM = p2mm.getCollisionDiameter(0);
        // MC = (MM + CC)/2
        return 2 * coreMC - coreMM;
    }

    public void boxLengthUpdated() {
        double L = box.getBoundary().getBoxSize().getX(1);
        double colloidSigma = getColloidSigma();
        p1WallColloid.setCollisionPosition(0, -L / 2 + colloidSigma / 2);
        p1WallColloid.setCollisionPosition(1, -L / 2 + rangeWallColloid);
        p1WallColloid.setCollisionPosition(2, +L / 2 - rangeWallColloid);
        p1WallColloid.setCollisionPosition(3, +L / 2 - colloidSigma / 2);

        double monomerSigma = p2mm.getCollisionDiameter(0);
        p1WallMonomer.setCollisionPosition(0, -L / 2 + monomerSigma / 2);
        p1WallMonomer.setCollisionPosition(1, -L / 2 + 2);
        p1WallMonomer.setCollisionPosition(2, +L / 2 - 2);
        p1WallMonomer.setCollisionPosition(3, +L / 2 - monomerSigma / 2);
    }

    public void setColloidSigma(double newColloidSigma) {
        if (newColloidSigma > 9) {
            throw new IllegalArgumentException("colloid sigma must be less than 9");
        }
        double L = box.getBoundary().getBoxSize().getX(1);
        p1WallColloid.setCollisionPosition(0, -L / 2 + newColloidSigma / 2);
        p1WallColloid.setCollisionPosition(3, +L / 2 - newColloidSigma / 2);

        p2mc.setCollisionDiameter(0, 0.5 * (p2mm.getCollisionDiameter(0) + newColloidSigma));
        configuration.setSigmaColloid(newColloidSigma);

        configuration.initializeCoordinates(box);
        if (integrator.getStepCount() > 0) {
            try {
                integrator.reset();
            } catch (ConfigurationOverlapException e) {
            }
        }
    }

}
