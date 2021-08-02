/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.ising;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPumpListener;
import etomica.graphics.*;
import etomica.integrator.IntegratorMCFasterer;
import etomica.lattice.LatticeCubicSimple;
import etomica.nbr.list.NeighborListManagerLattice;
import etomica.potential.BondingInfo;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.spin.heisenberg.P2Spin;


/**
 * Simulation of a simple 2D Ising model.  Prototype
 * for simulation of a more general magentic system.
 *
 * @author David Kofke
 */
public class IsingFasterer extends Simulation {

    private static final String APP_NAME = "Ising";
    public PotentialComputePairGeneral potentialMasterPair;
    public PotentialComputeField potentialMasterField;
    public PotentialComputeAggregate potentialMaster;
    public NeighborListManagerLattice nbrManager;
    public Box box;
    public SpeciesGeneral spins;
    public P2Spin potential;
    public P1MagneticFieldFasterer field;
    public MCMoveSpinFlipFasterer mcmove;
    public MeterSpinFasterer meter;
    public DataPumpListener pump;
    public AccumulatorAverageCollapsing dAcc;
    public IntegratorMCFasterer integrator;

    public IsingFasterer() {
        this(Space2D.getInstance(), 60);
    }

    /**
     *
     */
    public IsingFasterer(Space _space, int nCells) {
        super(_space);
        spins = SpeciesSpheresRotating.create(space, new ElementSimple(this));
        addSpecies(spins);

        box = this.makeBox();
        int numAtoms = space.powerD(nCells);
        box.setNMolecules(spins, numAtoms);
        box.getBoundary().setBoxSize(new Vector2D(nCells, nCells));
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(space, 1), space);
        config.initializeCoordinates(box);

        nbrManager = new NeighborListManagerLattice(getSpeciesManager(), box, 1, 1.1, BondingInfo.noBonding());
        nbrManager.setDoDownNeighbors(true);
        potentialMasterPair = new PotentialComputePairGeneral(getSpeciesManager(), box, nbrManager);
        potentialMasterField = new PotentialComputeField(getSpeciesManager(), box);
        potentialMaster = new PotentialComputeAggregate(potentialMasterPair, potentialMasterField);
        potential = new P2Spin(space);
        field = new P1MagneticFieldFasterer(space);
        integrator = new IntegratorMCFasterer(potentialMaster, random, 1.0, box);
        mcmove = new MCMoveSpinFlipFasterer(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcmove);

        getController().addActivity(new ActivityIntegrate(integrator));

        AtomType type = spins.getLeafType();
        potentialMasterField.setFieldPotential(type, field);
        potentialMasterPair.setPairPotential(type, type, potential);

        meter = new MeterSpinFasterer(box);
        dAcc = new AccumulatorAverageCollapsing();
        pump = new DataPumpListener(meter, dAcc, 10);
        dAcc.setPushInterval(10);
        integrator.getEventManager().addListener(pump);
    }

    public static void main(String[] args) {
        Space sp = Space2D.getInstance();
        IsingFasterer sim = new IsingFasterer(sp, 60);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
        ((SimulationRestart) simGraphic.getController().getReinitButton().getAction()).setConfiguration(null);
        IAction repaintAction = simGraphic.getPaintAction(sim.box);
        DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);

        simGraphic.remove(displayBox);
        displayBox.setBoxCanvas(new DisplayBoxSpin2DFasterer(displayBox, sp, sim.getController()));
        simGraphic.add(displayBox);
        DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integrator, "temperature");
        temperatureSlider.setMinimum(0.5);
        temperatureSlider.setMaximum(10.0);
        temperatureSlider.setShowBorder(true);
        simGraphic.add(temperatureSlider);
        temperatureSlider.setValue(sim.integrator.getTemperature());
        DeviceSlider fieldSlider = new DeviceSlider(sim.getController(), sim.field, "h");
        fieldSlider.setMinimum(-5.);
        fieldSlider.setMaximum(+5.);
        fieldSlider.setNMajor(5);
        fieldSlider.setValue(0.0);
        fieldSlider.setShowBorder(true);
        fieldSlider.setLabel("Magnetic field");
        simGraphic.add(fieldSlider);

        DisplayTextBoxesCAE boxes = new DisplayTextBoxesCAE();
        boxes.setAccumulator(sim.dAcc);
        boxes.setLabel("Magnetization");
        boxes.setLabelType(DisplayTextBox.LabelType.BORDER);
        simGraphic.add(boxes);

        simGraphic.getController().getReinitButton().setPostAction(repaintAction);
        simGraphic.makeAndDisplayFrame(APP_NAME);
    }
}
