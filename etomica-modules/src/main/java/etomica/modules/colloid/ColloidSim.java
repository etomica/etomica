/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.nbr.CriterionPositionWall;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.P2HardWrapper;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.random.RandomNumberGenerator;

/**
 * Colloid simulation.  Design by Alberto Striolo.
 *
 * @author Andrew Schultz
 */
public class ColloidSim extends Simulation {
    
    public PotentialMasterList potentialMaster;
    public SpeciesSpheresMono species, speciesColloid;
    public Box box;
    public IntegratorHard integrator;
    public P2HardWrapper potentialWrapper;
    public ActivityIntegrate activityIntegrate;
    public ConfigurationColloid configuration;
    public AtomLeafAgentManager<AtomArrayList> colloidMonomerBondManager;
    public AtomLeafAgentManager<AtomArrayList> monomerMonomerBondManager;
    public P2SquareWellMonomer p2mm;
    public P2HardSphereMC p2mc;
    public P2HardSphere p2pseudo;
    public int nGraft;
    public int chainLength;
    public P1Wall p1WallMonomer, p1WallColloid;
    public double epsWallWall = 0.25;
    public CriterionPositionWall criterionWallMonomer;
    
    public ColloidSim(Space _space) {
        super(_space);
        setRandom(new RandomNumberGenerator(1));

        //species
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        speciesColloid = new SpeciesSpheresMono(this, space);
        speciesColloid.setIsDynamic(true);
        addSpecies(speciesColloid);

        BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(null);
        potentialMaster = new PotentialMasterList(this, 6, boxAgentSource, new BoxAgentManager<NeighborCellManager>(boxAgentSource, this),
                new NeighborListManagerColloid.NeighborListAgentSourceColloid(6), _space);

        int nColloid = 1;
        chainLength = 50;
        nGraft = 12;

        double sigma = 1.0;
        double sigmaColloid = 7.5;
        double lambda = 1.5;
        double boxSize = 80;
        double epsMM = 1.0;

        //controller and integrator
        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTimeStep(0.02);
        integrator.setTemperature(2);
        integrator.setIsothermal(true);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);

        //potentials
        ((NeighborListManagerColloid) potentialMaster.getNeighborManager(box)).setSpeciesColloid(speciesColloid);
        ((NeighborListManagerColloid) potentialMaster.getNeighborManager(box)).setSpeciesMonomer(species);
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
        colloidMonomerBondManager = new AtomLeafAgentManager<AtomArrayList>(bondAgentSource, box);
        monomerMonomerBondManager = new AtomLeafAgentManager<AtomArrayList>(bondAgentSource, box);

        //instantiate several potentials for selection in combo-box
        p2mm = new P2SquareWellMonomer(space, monomerMonomerBondManager);
        p2mm.setCoreDiameter(sigma);
        p2mm.setLambda(lambda);
        p2mm.setBondFac(0.85);
        p2mm.setEpsilon(epsMM);
        potentialMaster.addPotential(p2mm, new AtomType[]{species.getLeafType(), species.getLeafType()});
        p2mc = new P2HardSphereMC(space, colloidMonomerBondManager);
        p2mc.setCollisionDiameter(0.5 * (sigma + sigmaColloid));
        p2mc.setBondFac(0.9);
        potentialMaster.addPotential(p2mc, new AtomType[]{species.getLeafType(), speciesColloid.getLeafType()});
        potentialMaster.setCriterion(p2mc, new CriterionNone() {
            public boolean needUpdate(IAtom atom) {
                return integrator.getStepCount() % 100 == 0;
            }
        });

        ((NeighborListManagerColloid) potentialMaster.getNeighborManager(box)).setPotentialMC(p2mc);

        p1WallMonomer = new P1Wall(space, monomerMonomerBondManager);
        p1WallMonomer.setBox(box);
        p1WallMonomer.setRange(2);
        p1WallMonomer.setSigma(1);
        p1WallMonomer.setEpsilon(Math.sqrt(epsMM * epsWallWall));
        potentialMaster.addPotential(p1WallMonomer, new AtomType[]{species.getLeafType()});
        criterionWallMonomer = new CriterionPositionWall(this);
        criterionWallMonomer.setBoundaryWall(true);
        criterionWallMonomer.setNeighborRange(3);
        criterionWallMonomer.setWallDim(1);
        potentialMaster.setCriterion(p1WallMonomer, criterionWallMonomer);

        p1WallColloid = new P1Wall(space, null);
        p1WallColloid.setBox(box);
        p1WallColloid.setRange(10);
        p1WallColloid.setSigma(7.5);
        p1WallColloid.setEpsilon(0.5);
        potentialMaster.addPotential(p1WallColloid, new AtomType[]{speciesColloid.getLeafType()});

        //construct box
        configuration = new ConfigurationColloid(space, species, speciesColloid, random);
        configuration.setNGraft(nGraft);
        configuration.setChainLength(chainLength);
        configuration.setSigmaColloid(sigmaColloid);
        configuration.setSigmaMonomer(sigma);
        configuration.setMonomerMonomerBondManager(monomerMonomerBondManager);
        configuration.setColloidMonomerBondManager(colloidMonomerBondManager);
        configuration.initializeCoordinates(box);

        p2pseudo = new P2HardSphere(space, 1, true);
        if (nGraft > 1) {
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
            p2pseudo.setCollisionDiameter(0.9 * Math.sqrt(minr2));
        }

        potentialMaster.addPotential(p2pseudo, new AtomType[]{species.getLeafType(), species.getLeafType()});
        potentialMaster.setCriterion(p2pseudo, new CriterionNone());

        ((NeighborListManagerColloid) potentialMaster.getNeighborManager(box)).setPotentialPseudo(p2pseudo);
        ((NeighborListManagerColloid) potentialMaster.getNeighborManager(box)).setChainLength(chainLength);

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        ColloidSim sim = new ColloidSim(space);
        sim.getController().actionPerformed();
    }
    
    public void setNumGraft(int newNumGraft) {
        if (nGraft == newNumGraft) return;
        configuration.setNGraft(newNumGraft);
        nGraft = newNumGraft;
        configuration.initializeCoordinates(box);

        if (newNumGraft > 1) {

            IAtom atom1 = box.getMoleculeList(species).get(0).getChildList().get(0);
            double minr2 = Double.POSITIVE_INFINITY;
            for (int j = chainLength; j<box.getMoleculeList(species).size(); j+=chainLength) {
                IAtom atom2 = box.getMoleculeList(species).get(j).getChildList().get(0);
                Vector dr = space.makeVector();
                dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                if (minr2 > r2) {
                    minr2 = r2;
                }
            }
            p2pseudo.setCollisionDiameter(0.9*Math.sqrt(minr2));
        }
        try {
            integrator.reset();
        }
        catch (ConfigurationOverlapException e) {}
    }

    public int getChainLength() {
        return chainLength;
    }
    
    public void setChainLength(int newChainLength) {
        if (chainLength == newChainLength) return;
        configuration.setChainLength(newChainLength);
        chainLength = newChainLength;
        configuration.initializeCoordinates(box);
        ((NeighborListManagerColloid)potentialMaster.getNeighborManager(box)).setChainLength(chainLength);

        try {
            integrator.reset();
        }
        catch (ConfigurationOverlapException e) {}
    }

    public double getColloidSigma() {
        return p1WallColloid.getSigma();
    }

    public void setColloidSigma(double newColloidSigma) {
        double oldSigma = p1WallColloid.getSigma();
        p1WallColloid.setRange(p1WallColloid.getRange()*newColloidSigma/oldSigma);
        p1WallColloid.setSigma(newColloidSigma);

        p2mc.setCollisionDiameter(0.5*(p2mm.getCoreDiameter()+newColloidSigma));
        configuration.setSigmaColloid(newColloidSigma);

        configuration.initializeCoordinates(box);
        if (integrator.getStepCount() > 0) {
            try {
                integrator.reset();
            }
            catch (ConfigurationOverlapException e) {}
        }
    }

    // reject everything.  we'll add them explicitly
    public static class CriterionNone implements NeighborCriterion {
        public boolean unsafe() {return false;}
        public void setBox(Box box) {}
        public void reset(IAtom atom) {}
        public boolean needUpdate(IAtom atom) {return false;}
        public boolean accept(IAtomList pair) {return false;}
    }
}
