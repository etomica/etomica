/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.PotentialMasterHybrid;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation for Adsorption module.
 * Design by Lev Gelb
 *
 * @author Andrew Schultz
 */
public class Adsorption extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono speciesA, speciesB;
    public final Box box;
    public final IntegratorHard integratorMD;
    public final IntegratorMC integratorMC;
    public final IntegratorHybrid integratorHybrid;
    public final ActivityIntegrate activityIntegrate;
    public final P2SquareWell p2AA, p2AB, p2BB;
    public final P1Wall p1WallA, p1WallB;
    public final MyMCMove mcMoveIDA, mcMoveIDB;
    
    public Adsorption(Space _space) {
        super(_space);
        //construct box
        box = new Box(new BoundaryRectangularSlit(1, 20.0, space), space);
        addBox(box);
        PotentialMasterHybrid potentialMaster = new PotentialMasterHybrid(this, 2, space); //List(this, 2.0);

        //controller and integrator
        integratorMD = new IntegratorHard(this, potentialMaster, space, box);
        integratorMD.setTimeStep(0.005);
        integratorMD.setIsothermal(true);
        integratorMD.setThermostatInterval(10000);

        integratorMC = new IntegratorMC(potentialMaster, random, 2, box);
        integratorMC.setTemperature(1);

        integratorHybrid = new IntegratorHybrid(potentialMaster, integratorMD, integratorMC, 2);

        activityIntegrate = new ActivityIntegrate(integratorHybrid);
        getController().addAction(activityIntegrate);

        double sigma = 1;
        double lambda = 1.5;
        double epsilon = 1.0;
        double epsilonWF = 5.0;

        //species and potentials
        speciesA = new SpeciesSpheresMono(space, new ElementSimple(this));
        speciesA.setIsDynamic(true);
        addSpecies(speciesA);
        speciesB = new SpeciesSpheresMono(space, new ElementSimple(this));
        speciesB.setIsDynamic(true);
        addSpecies(speciesB);


        mcMoveIDA = new MyMCMove(integratorMC, random, space, 0.9, sigma, 1);
        mcMoveIDA.setMu(-12);
        integratorMC.getMoveManager().addMCMove(mcMoveIDA);
        mcMoveIDA.setSpecies(speciesA);
        mcMoveIDA.setBox(box);

        mcMoveIDB = new MyMCMove(integratorMC, random, space, 0.9, sigma, 1);
        mcMoveIDB.setMu(-Double.POSITIVE_INFINITY);
        integratorHybrid.setMCMoveInsertDelete(mcMoveIDA, mcMoveIDB);
        mcMoveIDB.setSpecies(speciesB);
        mcMoveIDB.setBox(box);


        p2AA = new P2SquareWell(space, sigma, lambda, epsilon, false);
        potentialMaster.addPotential(p2AA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
        p2AB = new P2SquareWell(space, sigma, lambda, epsilon, false);
        potentialMaster.addPotential(p2AB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
        p2BB = new P2SquareWell(space, sigma, lambda, epsilon, false);
        potentialMaster.addPotential(p2BB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});


        p1WallA = new P1Wall(space);
        p1WallA.setSigma(sigma);
        p1WallA.setRange(sigma / 2);
        p1WallA.setEpsilon(epsilonWF);
        p1WallA.setThermalize(integratorMC, 0.0, random);

        p1WallB = new P1Wall(space);
        p1WallB.setSigma(sigma);
        p1WallB.setRange(sigma / 2);
        p1WallB.setEpsilon(epsilonWF);
        p1WallB.setThermalize(integratorMC, 0.0, random);

        potentialMaster.addPotential(p1WallA, new AtomType[]{speciesA.getLeafType()});
        potentialMaster.addPotential(p1WallB, new AtomType[]{speciesB.getLeafType()});

        integratorMD.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        Vector dim = space.makeVector();
        dim.E(8 * sigma);
        dim.setX(1, 12 * sigma);
        box.getBoundary().setBoxSize(dim);

        box.setNMolecules(speciesA, 40);

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        
        Adsorption sim = new Adsorption(space);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "Catalysis", 1);
        simGraphic.makeAndDisplayFrame();
    }
}
