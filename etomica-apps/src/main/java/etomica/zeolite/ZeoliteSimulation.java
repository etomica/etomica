/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;

import etomica.action.SimulationRestart;

import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.*;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2WCA;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;


/**
 * 
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.  
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class ZeoliteSimulation extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    private static final long serialVersionUID = 1L;

    private static final String APP_NAME = "Zeolite Simulation";

    /**
     * The Box holding the atoms. 
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    //public final IntegratorHard integrator;
    public final IntegratorVelocityVerlet integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono[] species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    //public final P2HardSphere potential;
    public final P2LennardJones potentialMM;

    public DisplayPlot ePlot;
    private int nAtomsMeth;
    private int interval;
    private SpeciesSpheresMono sp;
    private String filename;

    public ZeoliteSimulation() {

        super(Space3D.getInstance());

        //Additions for Zeolite Calculations
        //Start by reading the first line, which is number of Atoms
        //String fileName = "pbu2";
        ConfigurationFileXYZ config = new ConfigurationFileXYZ("2unitcell.xyz", space);
        int[] numAtoms = config.getNumAtoms();

        species = new SpeciesSpheresMono[numAtoms.length];
        for (int i = 0; i < numAtoms.length; i++) {
            species[i] = new SpeciesSpheresMono(this, space);
            species[i].setIsDynamic(true);
            addSpecies(species[i]);
            if (i != (numAtoms.length - 1)) {
                // all elements except the last (methane) are fixed
                ((ElementSimple) (species[i].getLeafType()).getElement()).setMass(Double.POSITIVE_INFINITY);
            } else {
                ((ElementSimple) species[i].getLeafType().getElement()).setMass(16);
            }
        }

        nAtomsMeth = numAtoms[numAtoms.length - 1];
        double neighborRangeFac = 1.2;
        //defaults.makeLJDefaults();
        //Setting sizes of molecules

        PotentialMasterList potentialMaster = new PotentialMasterList(this, 1.6, space);
        double range = 8.035;
        potentialMaster.setRange(3.214 * neighborRangeFac * 2.5);


        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(10);
        integrator.setTimeStep(0.00611);
        integrator.setTemperature(Kelvin.UNIT.toSim(298.0));


        getController2().addActivity(new ActivityIntegrate2(integrator, 2, true), 1000000, 2);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        for (int i = 0; i < numAtoms.length; i++) {
            box.setNMolecules(species[i], numAtoms[i]);
        }
        //Setting up potential for Methane-Methane interactions
        potentialMM = new P2LennardJones(space);
        //Setting sigma and epsilon (from J. Chem. Soc Faraday Trans 1991
        potentialMM.setEpsilon(Kelvin.UNIT.toSim(147.95));
        //potentialMM.setEpsilon(1);
        potentialMM.setSigma(3.0);
        //Setting up potential for Methane-Oxygen interactions
        P2LennardJones potentialMO = new P2LennardJones(space);
        //Setting sigma and epsilon
        potentialMO.setEpsilon(Kelvin.UNIT.toSim(133.3));
        //potentialMO.setEpsilon(100);
        potentialMO.setSigma(3.214);

        //Setting up Methane - Silicon interactions
        //P2LennardJones potentialMS = potentialMO;
        P2WCA potentialMS = new P2WCA(space, 1.18, potentialMO.getEpsilon());

        //Wrap LJ potentials to truncate
        P2SoftSphericalTruncated MM = new P2SoftSphericalTruncated(space, potentialMM, 2.5 * potentialMM.getSigma());
        P2SoftSphericalTruncated MO = new P2SoftSphericalTruncated(space, potentialMO, 2.5 * potentialMO.getSigma());
        //P2SoftSphericalTruncated MS = new P2SoftSphericalTruncated(potentialMS,2.5*potentialMS.getSigma());


        potentialMaster.addPotential(MM, new AtomType[]{species[2].getLeafType(), species[2].getLeafType()});
        potentialMaster.addPotential(MO, new AtomType[]{species[0].getLeafType(), species[2].getLeafType()});
        potentialMaster.addPotential(potentialMS, new AtomType[]{species[1].getLeafType(), species[2].getLeafType()});

        //Initializes the coordinates and positions
        config.initializeCoordinates(box);
        box.getBoundary().setBoxSize(config.getUpdatedDimensions());
        //integrator.addListener(new BoxImposePbc(box));

        //PARAMETERS For Simulation Run
        //activityIntegrate.setMaxSteps(5000000);
        double ts = 0.00611;
        integrator.setTimeStep(ts);
        interval = 2000;
        integrator.setThermostatInterval(interval / 1000);

        //      Adding coordinate writer by Mike Sellars

        filename = (numAtoms[2] + "_" + getController2().getMaxSteps() + "_" + ts + "_" + interval + "_WCA");
        sp = species[2];
        /*
        MSDCoordWriter coordWriter = new MSDCoordWriter(this.space, filename,sp);
        coordWriter.setBox(this.box);
        coordWriter.setNatoms(numAtoms[2]);
        coordWriter.setIntegrator(integrator);
        coordWriter.setWriteInterval(interval);
        */
    } //end of constructor

    public static void Converter(String inputFile) {
		// TODO Auto-generated method stub
		String outputFile = inputFile+"__Result.txt";
		MSDProcessor proc = new MSDProcessor(Space3D.getInstance(),inputFile,outputFile);

        //proc.setDeltaTmax(1);
		proc.fillArrays();
		System.out.println("Converter done");
	}

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        //defaults.temperature = Kelvin.UNIT.toSim(298.0);
        ZeoliteSimulation sim = new ZeoliteSimulation();
        zeoliteSimGraphic simGraphic = new zeoliteSimGraphic(sim, sim.space, APP_NAME);
        int num = sim.species.length;
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpecies(sim.species[num-1]);
        nSelector.setBox(sim.box);
        simGraphic.add(nSelector);

        //Energy
        int history = sim.getInterval()*10;
        //Settings

        MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.getHistory().setHistoryLength(history);
        AccumulatorAverageCollapsing enAcc = new AccumulatorAverageCollapsing();
        enAcc.setPushInterval(20);
        DataFork enFork = new DataFork(new IDataSink[]{energyHistory, enAcc});
        DataPump energyPump = new DataPump(eMeter, enFork);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(energyPump);
        pumpListener.setInterval(10);
        sim.integrator.getEventManager().addListener(pumpListener);
        energyHistory.setPushInterval(10);
        simGraphic.getController().getDataStreamPumps().add(energyPump);

        /*
        MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.potentialMaster);
        peMeter.setBox(sim.box);
        AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setHistoryLength(history);
        AccumulatorAverage peAccumulator = new AccumulatorAverage(sim);
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        IntervalActionAdapter peAdapter = new IntervalActionAdapter(pePump);
        peAdapter.setActionInterval(10);
        peHistory.setPushInterval(10);
        sim.register(peMeter,pePump);
        sim.integrator.addListener(peAdapter);

		MeterKineticEnergy keMeter = new MeterKineticEnergy();
        keMeter.setBox(sim.box);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setHistoryLength(history);
        DataPump kePump = new DataPump(keMeter, keHistory);
        IntervalActionAdapter keAdapter = new IntervalActionAdapter(kePump);
        keAdapter.setActionInterval(10);
        keHistory.setPushInterval(10);
        sim.register(keMeter,kePump);
        sim.integrator.addListener(keAdapter);


        MeterTemperature temp = new MeterTemperature();
		temp.setBox(sim.box);
		AccumulatorHistory tHistory = new AccumulatorHistory();
        tHistory.setHistoryLength(history);
        AccumulatorAverage tAccumulator = new AccumulatorAverage(sim);
        DataFork tFork = new DataFork(new DataSink[]{tHistory, tAccumulator});
        DataPump tPump = new DataPump(temp, tFork);
        IntervalActionAdapter tAdapter = new IntervalActionAdapter(tPump);
        tAdapter.setActionInterval(10);
        tHistory.setPushInterval(10);
        sim.register(temp,tPump);
        sim.integrator.addListener(tAdapter);
        */

        DisplayPlot ePlot = sim.ePlot;
        ePlot = new DisplayPlot();
        ePlot.setLabel("Energy");
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        //peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        //keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
		//tHistory.setDataSink(ePlot.getDataSet().makeDataSink());
		ePlot.setDoLegend(true);
		simGraphic.add(ePlot);


        simGraphic.makeAndDisplayFrame(APP_NAME);

        ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box).getColorScheme();
        DiameterHashByType diameterManager = (DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash();
        double[] atomicSize = new double[3];
        atomicSize[0] = 0.73;
        atomicSize[1] = 1.18;
        atomicSize[2] = 2.088;
        for(int i=0;i<sim.species.length;i++){
        	switch(i){
        		case 0:
        			colorScheme.setColor(sim.species[i].getLeafType(), java.awt.Color.red);
        			break;
        		case 1:
        			colorScheme.setColor(sim.species[i].getLeafType(), java.awt.Color.blue);
        			break;
        		default:
        			colorScheme.setColor(sim.species[i].getLeafType(), java.awt.Color.white);
        	}

            diameterManager.setDiameter(sim.species[i].getLeafType(), atomicSize[i]);
        }

    }//end of main

    int getInterval() {
        return interval;
    }

    int getMethane() {
        return nAtomsMeth;
    }

    String getFileName() {
        return filename;
    }

    SpeciesSpheresMono getSpeciesRMS() {
        return sp;
    }

}//end of class
