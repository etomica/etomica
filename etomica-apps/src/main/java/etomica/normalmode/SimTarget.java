/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

import java.io.FileWriter;
import java.io.IOException;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * @author Andrew Schultz
 */
public class SimTarget extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorHard integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public BravaisLattice lattice;
    public Primitive primitive;
    public CoordinateDefinition coordinateDefinition;
    public SimTarget(Space _space, int numAtoms, double density) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        PotentialMaster potentialMaster = (space.D() == 1 ? new PotentialMasterList(this, space) : new PotentialMasterMonatomic(this));

        Potential potential = new P2HardSphere(space, 1.0, false);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType,
                sphereType});
        int nCells;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0 / density);
            boundary = new BoundaryRectangularPeriodic(space, numAtoms / density);
            integrator.setNullPotential(new P1HardPeriodic(space), sphereType);
            nCells = numAtoms;
        } else {
            primitive = new PrimitiveFcc(space, 1);
            double v = primitive.unitCell().getVolume();
            primitive.scaleSize(Math.pow(v * density, -1.0 / 3.0));
            nCells = (int) Math.round(Math.pow(numAtoms, 1.0 / 3.0));
            boundary = new BoundaryDeformableLattice(primitive, new int[]{nCells, nCells, nCells});
        }
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorHard(this, potentialMaster, box);

        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(integrator);
        double timeStep = 0.4;
        integrator.setTimeStep(timeStep);
        getController().addAction(activityIntegrate);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, space);
        coordinateDefinition.initializeCoordinates(new int[]{nCells, nCells, nCells});

        if (potentialMaster instanceof PotentialMasterList) {
            double neighborRange;
            if (space.D() == 1) {
                neighborRange = 1.01 / density;
            } else {
                //FCC
                double L = Math.pow(4.01 / density, 1.0 / 3.0);
                neighborRange = L / Math.sqrt(2.0);
            }
            ((PotentialMasterList) potentialMaster).setRange(neighborRange);
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMaster).getNeighborManager(box).reset();
        }
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        //set up simulation parameters
        int D = 3;
        int numMolecules = 27;
        double density = 1.30;
        double harmonicFudge = .1;
        double simTime = 1000;
        double temperature = 1;
        if (D == 1) {
            numMolecules = 3;
            density = 0.5;
            simTime = 1000000;
        }
        String filename = "normal_modes"+D+"D_"+numMolecules+"_130";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            numMolecules = Integer.parseInt(args[3]);
        }
        if (args.length > 4) {
            harmonicFudge = Double.parseDouble(args[4]);
        }


        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" hard sphere simulation, measuring harmonic energy");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println(simTime+" time units");
        System.out.println("output data to "+filename);

        //instantiate simulation
        SimTarget sim = new SimTarget(Space.getInstance(D), numMolecules, density);

        NormalModes normalModes = null;
        if(D == 1) {
            normalModes = new NormalModes1DHR(sim.boundary, numMolecules);
        } else {
            normalModes = new NormalModesFromFile(filename, D);
        }
        normalModes.setHarmonicFudge(harmonicFudge);
        normalModes.setTemperature(temperature);

        //add meters to for FEP averages
        //this one does averaging of total energy and its Boltzmann factor
        //so long as we're using a MeterHarmonicSingleEnergy, we'll use that instead to get the sum
//        MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy(new CoordinateDefinitionLeaf(sim.getSpace()), normalModes);
//        harmonicEnergy.setBox(sim.box);
//        DataPump foo = new DataPump(harmonicEnergy, null);
//        IntervalActionAdapter bar = new IntervalActionAdapter(foo, sim.integrator);
//        bar.setActionInterval(50);

        //this one does averaging of Boltzmann factors of each mode
        MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy(sim.coordinateDefinition, normalModes);
        harmonicSingleEnergy.setBox(sim.box);
        DataPump pump = new DataPump(harmonicSingleEnergy, null);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(50);

        DataFork harmonicSingleFork = new DataFork();
        pump.setDataSink(harmonicSingleFork);
        BoltzmannProcessor boltz = new BoltzmannProcessor();
        boltz.setTemperature(temperature);
        harmonicSingleFork.addDataSink(boltz);
        DataFork boltzFork = new DataFork();
        boltz.setDataSink(boltzFork);
//        DataProcessorCorrelationMatrix boltzCorrelation = new DataProcessorCorrelationMatrix();
//        boltzFork.addDataSink(boltzCorrelation);
        AccumulatorAverage harmonicSingleAvg = new AccumulatorAverageFixed(10);
        boltzFork.addDataSink(harmonicSingleAvg);

        DataProcessorSum summer = new DataProcessorSum();
        harmonicSingleFork.addDataSink(summer);
        DataFork harmonicFork = new DataFork();
        summer.setDataSink(harmonicFork);
        AccumulatorAverage harmonicAvg = new AccumulatorAverageFixed(10);
        harmonicFork.addDataSink(harmonicAvg);
        boltz = new BoltzmannProcessor();
        boltz.setTemperature(temperature);
        harmonicFork.addDataSink(boltz);
        AccumulatorAverage harmonicBoltzAvg = new AccumulatorAverageFixed(10);
        boltz.setDataSink(harmonicBoltzAvg);

        //start simulation
        int nSteps = (int) (simTime / sim.integrator.getTimeStep());
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();

        //get averages and confidence limits for harmonic energy
        double avgHarmonicEnergy = ((DataDouble) ((DataGroup) harmonicAvg.getData()).getData(harmonicAvg.AVERAGE.index)).x;
        double errorHarmonicEnergy = ((DataDouble) ((DataGroup) harmonicAvg.getData()).getData(harmonicAvg.ERROR.index)).x;
        System.out.println("avg harmonic energy: "+avgHarmonicEnergy+" +/- "+errorHarmonicEnergy);

        //compute free-energy quantities, independent-mode approximation
        DataDoubleArray harmonicModesAvg = (DataDoubleArray) ((DataGroup) harmonicSingleAvg.getData()).getData(harmonicSingleAvg.AVERAGE.index);
        DataDoubleArray harmonicModesErr = (DataDoubleArray) ((DataGroup) harmonicSingleAvg.getData()).getData(harmonicSingleAvg.ERROR.index);
        double deltaA = 0;
        double deltaAerr = 0;
        int nData = harmonicModesAvg.getLength();
        for (int i=0; i<nData; i++) {
            deltaA += Math.log(harmonicModesAvg.getValue(i));
            deltaAerr += harmonicModesErr.getValue(i)/harmonicModesAvg.getValue(i);
        }
        System.out.println("Harmonic free energy correction (independent approx): "+deltaA+" +/- "+deltaAerr);

        double[][] omega2 = normalModes.getOmegaSquared();
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        double AHarmonic = 0;
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));//coeffs in log?
            }
        }
        int totalCells = numMolecules;
        int basisSize = 1;
        if (totalCells % 2 == 0) {
            AHarmonic += Math.log(Math.pow(2.0, basisSize*D*(totalCells - Math.pow(2,D))/2.0) / Math.pow(totalCells,0.5*D));
        }
        else {
            AHarmonic += Math.log(Math.pow(2.0, basisSize*D*(totalCells - 1)/2.0) / Math.pow(totalCells,0.5*D));
        }

        //results for averaging without independent-mode approximation
        deltaA = ((DataDouble) ((DataGroup) harmonicBoltzAvg.getData()).getData(harmonicBoltzAvg.AVERAGE.index)).x;
        deltaAerr = ((DataDouble) ((DataGroup) harmonicBoltzAvg.getData()).getData(harmonicBoltzAvg.ERROR.index)).x / deltaA;
        deltaA = Math.log(deltaA);

        System.out.println("Harmonic free energy correction: "+deltaA+" +/- "+deltaAerr);
        System.out.println("Harmonic free energy correction per atom: "+deltaA/numMolecules+" +/- "+deltaAerr/numMolecules);

        System.out.println("Harmonic-reference free energy: "+AHarmonic);

        if(D==1) {
            double AHR = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules) + SpecialFunctions.lnFactorial(numMolecules) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }

        try {
            // write averages of exp(-u/kT) for each normal mode
            FileWriter fileWriterE = new FileWriter(filename+".e");
            int[] idx = new int[2];
            for (int i=0; i<harmonicModesAvg.getArrayShape(0); i++) {
                idx[0] = i;
                idx[1] = 0;
                fileWriterE.write(Double.toString(harmonicModesAvg.getValue(idx)));
                for (int j=1; j<harmonicModesAvg.getArrayShape(1); j++) {
                    idx[1] = j;
                    fileWriterE.write(" "+harmonicModesAvg.getValue(idx));
                }
                fileWriterE.write("\n");
            }
            fileWriterE.close();
        }
        catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data "+e);
        }

//        DataDoubleArray matrix = (DataDoubleArray)boltzCorrelation.getData();
//        try {
//            FileWriter fileWriterEE = new FileWriter(filename+".ee");
//            int k=0;
//            for (int i=0; i<matrix.getArrayShape(0); i++) {
//                for (int j=0; j<matrix.getArrayShape(1); j++) {
//                    fileWriterEE.write(i+" "+j+" "+(matrix.getValue(k)/(harmonicModesAvg.getValue(i)*harmonicModesAvg.getValue(j))-1)+"\n");
//                    k++;
//                }
//            }
//            fileWriterEE.close();
//        }
//        catch (IOException e) {
//            throw new RuntimeException("Oops, failed to write data "+e);
//        }
//        nData = matrix.getLength();
//        deltaA = 0;
//        for (int i=0; i<nData; i++) {
//            deltaA += Math.log(matrix.getValue(i));
//        }
//        System.out.println("pair approximation to free energy correction "+deltaA/(numMolecules-1));
    }
}
