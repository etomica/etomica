/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.*;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Energy;
import etomica.units.SimpleUnit;
import etomica.util.ParameterBase;

/**
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 */
public class HSMD3D extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    private static final long serialVersionUID = 1L;
    /**
     * The Box holding the atoms.
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardSphere potential;

    public final PotentialMaster potentialMaster;

    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public HSMD3D(Space _space) {
        this(_space, new HSMD3DParam());
    }

    public HSMD3D(Space _space, HSMD3DParam params) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(_space);

        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, 3.0, space) : new PotentialMasterMonatomic(this);

        int numAtoms = params.nAtoms;
        double neighborRangeFac = 1.5;
        double sigma = 1.0;
        if (params.useNeighborLists) {
            ((PotentialMasterList) potentialMaster).setRange(neighborRangeFac * sigma);
        }

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);
        integrator.setTemperature(2.0);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        potential = new P2HardSphere(space, sigma, false);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.eta * 2 * space.D() / Math.PI);
        inflater.actionPerformed();
        if (space.D() == 3) {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        } else {
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        }
        //deformed
//        box.setBoundary(
//            new etomica.space.BoundaryDeformablePeriodic(
//            space,getRandom(),
//            new IVector[]{
//              new Vector3D(-4,1,1),
//              new Vector3D(2,6,4),
//              new Vector3D(1,2,6)}));
        //truncated octahedron
//        box.setBoundary(
//            new etomica.space3d.BoundaryTruncatedOctahedron(this));

        integrator.setBox(box);

        if (params.useNeighborLists) {
            NeighborListManager nbrManager = ((PotentialMasterList) potentialMaster).getNeighborManager(box);
            integrator.getEventManager().addListener(nbrManager);
        } else {
            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        }
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final String APP_NAME = "HSMD3D";

        Space sp = Space3D.getInstance();
        HSMD3DParam params = new HSMD3DParam();
        final HSMD3D sim = new HSMD3D(sp, params);
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, sim.space, sim.getController());
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim, sp, sim.getController()));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);

        nSelector.setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(nSelector);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);

        PotentialMasterCell potentialMasterSW = new PotentialMasterCell(sim, 2.0, sim.getSpace());
        potentialMasterSW.setCellRange(2);
        potentialMasterSW.getNbrCellManager(sim.box).setDoApplyPBC(true);

        double lambda = 2;

        P2SquareWell potentialSW = new etomica.potential.P2SquareWell(sim.getSpace(), 1.0, lambda, 1.0, false);

        potentialMasterSW.addPotential(potentialSW, new AtomType[]{sim.species.getLeafType(), sim.species.getLeafType()});
        MeterPotentialEnergy meterPESW = new MeterPotentialEnergy(potentialMasterSW) {
            public double getDataAsScalar() {
                ((PotentialMasterCell) potential).getNbrCellManager(box).assignCellAll();
                return super.getDataAsScalar();
            }
        };
        meterPESW.setBox(sim.box);

        AccumulatorHistory uSWHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        DataPumpListener pumpSW = new DataPumpListener(meterPESW, uSWHistory, 10);
        sim.integrator.getEventManager().addListener(pumpSW);
        simGraphic.getController().getDataStreamPumps().add(pumpSW);
        DisplayPlot plotPE = new DisplayPlot();
        uSWHistory.setDataSink(plotPE.getDataSet().makeDataSink());
        plotPE.setLabel("SW");
        plotPE.setUnit(new SimpleUnit(Energy.DIMENSION, params.nAtoms, "", "", false));
        simGraphic.add(plotPE);
    }

    public static HSMD3DParam getParameters() {
        return new HSMD3DParam();
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class HSMD3DParam extends ParameterBase {
        public int nAtoms = 12800;
        public double eta = 0.45;
        public boolean useNeighborLists = true;
    }
}
