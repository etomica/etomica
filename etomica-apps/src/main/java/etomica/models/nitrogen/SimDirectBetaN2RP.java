/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.ActivityIntegrate2;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * Direct Sampling for Rotational Perturbation
 *   
 * @author Tai Boon Tan
 */
public class SimDirectBetaN2RP extends Simulation {

    public SimDirectBetaN2RP(Space space, int numMolecules, double density, double temperature, double[] angle) {
        super(space);

        species = new SpeciesN2(space);
        addSpecies(species);

        PotentialMaster potentialMasterTarg = new PotentialMaster();
        PotentialMaster potentialMasterRef = new PotentialMaster();

        // TARGET
        double ratio = 1.631;
        double aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double cDim = aDim * ratio;
        System.out.println("\naDim: " + aDim + " ;cDim: " + cDim);
        int nC = (int) Math.pow(numMolecules / 1.999999999, 1.0 / 3.0);

        Basis basisHCP = new BasisHcp();
        BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nC, nC, nC});

        Vector[] boxDim = new Vector[3];
        boxDim[0] = Vector.of(nC * aDim, 0, 0);
        boxDim[1] = Vector.of(-nC * aDim * Math.cos(Degree.UNIT.toSim(60)), nC * aDim * Math.sin(Degree.UNIT.toSim(60)), 0);
        boxDim[2] = Vector.of(0, 0, nC * cDim);
        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        boxTarg = this.makeBox(boundary);
        boxTarg.setNMolecules(species, numMolecules);


        int[] nCells = new int[]{1, 1, 1};
        Primitive primitive = new PrimitiveHexagonal(space, nC * aDim, nC * cDim);

        coordinateDefTarg = new CoordinateDefinitionNitrogen(this, boxTarg, primitive, basis, space);
        coordinateDefTarg.setIsBeta();
        coordinateDefTarg.setOrientationVectorBeta(space);
        coordinateDefTarg.initializeCoordinates(nCells);

        double rCScale = 0.475;
        double rc = aDim * nC * rCScale;
        System.out.println("Truncation Radius (" + rCScale + " Box Length): " + rc);
        P2Nitrogen potentialTarg = new P2Nitrogen(space, rc);
        potentialTarg.setBox(boxTarg);

        PRotConstraint pRotConstraintTarg = new PRotConstraint(space, coordinateDefTarg, boxTarg);
        pRotConstraintTarg.setConstraintAngle(angle[0]);

        potentialMasterTarg.addPotential(potentialTarg, new ISpecies[]{species, species});
        potentialMasterTarg.addPotential(pRotConstraintTarg, new ISpecies[]{species});

        MCMoveMoleculeCoupled moveTarg = new MCMoveMoleculeCoupled(potentialMasterTarg, getRandom(), space);
        moveTarg.setBox(boxTarg);
        moveTarg.setPotential(potentialTarg);

        MCMoveRotateMolecule3D rotateTarg = new MCMoveRotateMolecule3D(potentialMasterTarg, getRandom(), space);
        rotateTarg.setBox(boxTarg);

        integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature, boxTarg);
        integratorTarg.getMoveManager().addMCMove(moveTarg);
        integratorTarg.getMoveManager().addMCMove(rotateTarg);

        MeterPotentialEnergy meterPETarg = new MeterPotentialEnergy(potentialMasterTarg, boxTarg);
        System.out.println("lattice energy (sim unit): " + meterPETarg.getDataAsScalar());

        // Reference System

        PRotConstraint pRotConstraintRef = new PRotConstraint(space, coordinateDefTarg, boxTarg);
        pRotConstraintRef.setConstraintAngle(angle[1]);

        potentialMasterRef.addPotential(potentialTarg, new ISpecies[]{species, species});
        potentialMasterRef.addPotential(pRotConstraintRef, new ISpecies[]{species});

        MeterPotentialEnergy meterPERef = new MeterPotentialEnergy(potentialMasterRef, boxTarg);

        MeterBoltzmannDirect meterBoltzmann = new MeterBoltzmannDirect(integratorTarg, meterPERef);
        boltzmannAverage = new AccumulatorAverageFixed(100);

        DataPump boltzmannPump = new DataPump(meterBoltzmann, boltzmannAverage);
        IntegratorListenerAction boltzmannPumpListener = new IntegratorListenerAction(boltzmannPump, 100);
        integratorTarg.getEventManager().addListener(boltzmannPumpListener);

        activityIntegrate = new ActivityIntegrate(integratorTarg);
        getController().addAction(activityIntegrate);
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimDirectBetaN2RP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        
        double density = params.density;
        double[] angle = params.angle;
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
  
        System.out.println("Running beta-phase Nitrogen RP direct sampling simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.println("perturbing from angle=" + angle[0] + " into " +angle[1]);
        System.out.println("with simulation steps of " + numSteps);
        
        SimDirectBetaN2RP sim = new SimDirectBetaN2RP(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), angle);

        MeterOrientationDistribution meterOrient = new MeterOrientationDistribution(sim.boxTarg, sim.coordinateDefTarg, sim.species);
        IntegratorListenerAction meterOrientListener = new IntegratorListenerAction(meterOrient);
        meterOrientListener.setInterval(numMolecules);                                      
        sim.integratorTarg.getEventManager().addListener(meterOrientListener);    
        
        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorTarg), numSteps/10);
System.out.println("equilibration finished");
        sim.getController().reset();


        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorTarg), numSteps);

        double average = sim.boltzmannAverage.getData().getValue(sim.boltzmannAverage.AVERAGE.index);
        double error = sim.boltzmannAverage.getData().getValue(sim.boltzmannAverage.ERROR.index);
        meterOrient.writeUdistribution("rotDistA"+angle[0]);
        
        System.out.println("boltzmann average "+angle[0] +" to "+angle[1]+": " + average + " ;err: " + error);
        
        long endTime = System.currentTimeMillis();
        System.out.println("\nEnd Time: " + endTime);
        System.out.println("Time taken (s): " + (endTime - startTime)/1000);
       
    }

    private static final long serialVersionUID = 1L;
    protected ActivityIntegrate activityIntegrate;
    protected AccumulatorAverageFixed boltzmannAverage;
    protected Box boxTarg;
    protected SpeciesN2 species;
    protected CoordinateDefinitionNitrogen coordinateDefTarg;
    protected IntegratorMC integratorTarg;
    
    /**
     * Inner class for parameters understood by the SimOverlapBetaN2RP constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 128;
        public double density = 0.025;
        public double[] angle = new double[]{86, 85.5};
        public int D = 3;
        public long numSteps =500000;
        public double temperature = 40;
    }
}
