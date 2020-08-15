/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.ising;

import etomica.action.IAction;
import etomica.action.SimulationRestart;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.nbr.site.NeighborSiteManager;
import etomica.nbr.site.PotentialMasterSite;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;


/**
 * Simulation of a simple 2D Ising model.  Prototype
 * for simulation of a more general magentic system.
 *
 * @author David Kofke
 */
public class Ising extends Simulation {

    private static final String APP_NAME = "Ising";
    private static final long serialVersionUID = 2L;
    public PotentialMasterSite potentialMaster;
    public Box box;
    public SpeciesSpheresMono spins;
    public P2Spin potential;
    public P1MagneticField field;
    public MCMoveSpinFlip mcmove;
    public MeterSpin meter;
    public DataPump pump;
    public AccumulatorAverageCollapsing dAcc;
    private IntegratorMC integrator;

    public Ising() {
        this(Space2D.getInstance(), 60);
    }

    /**
     *
     */
    public Ising(Space _space, int nCells) {
        super(_space);
        spins = new SpeciesSpheresMono(this, space);
        addSpecies(spins);
        potentialMaster = new PotentialMasterSite(this, nCells, space);
        box = this.makeBox();
        int numAtoms = space.powerD(nCells);
        box.setNMolecules(spins, numAtoms);
        new ConfigurationAligned().initializeCoordinates(box);

        potential = new P2Spin(space);
        field = new P1MagneticField(space);
        integrator = new IntegratorMC(this, potentialMaster, box);
        mcmove = new MCMoveSpinFlip(potentialMaster, getRandom());
        integrator.getMoveManager().addMCMove(mcmove);

        getController().addActivity(new ActivityIntegrate(integrator));

        AtomType type = spins.getLeafType();
        potentialMaster.addPotential(field, new AtomType[]{type});
        potentialMaster.addPotential(potential, new AtomType[]{type, type});

        meter = new MeterSpin(space);
        meter.setBox(box);
        dAcc = new AccumulatorAverageCollapsing();
        pump = new DataPump(meter, dAcc);
        dAcc.setPushInterval(10);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(10);
        integrator.getEventManager().addListener(pumpListener);
    }

    public static void main(String[] args) {
        Space sp = Space2D.getInstance();
        Ising sim = new Ising(sp, 60);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
        ((SimulationRestart) simGraphic.getController().getReinitButton().getAction()).setConfiguration(null);
        IAction repaintAction = simGraphic.getPaintAction(sim.box);
        DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);

        simGraphic.remove(displayBox);
        NeighborSiteManager neighborSiteManager = (NeighborSiteManager)sim.potentialMaster.getBoxCellManager(sim.box);
        displayBox.setBoxCanvas(new DisplayBoxSpin2D(displayBox,neighborSiteManager, sp, sim.getController()));
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
