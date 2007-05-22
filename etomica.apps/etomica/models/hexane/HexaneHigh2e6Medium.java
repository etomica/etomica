package etomica.models.hexane;

import etomica.integrator.IntervalActionAdapter;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
/**
 * @author nancycribbin
 *  
 */

/*
 * We use a PotentialMaster, rather than a PotentialMasterNbr, so that we do not
 * need to deal with cells, which BoundaryDeformablePeriodic cannot deal with at
 * this time.
 * 
 * @author nancycribbin
 * 
 */

public class HexaneHigh2e6Medium extends Simulation {

    private static final long serialVersionUID = 1L;

    public static void main(String[] args) {
        //defaults
        int D = 3;
        int nA = 36;
        double density = 0.40;
        long nSteps = 2000000;
        int numMolecules = nA; //144
        
        String filename = "normal_modes"+nA+"_"+((int)(density*100))+"hexane"+nSteps;

        System.out.println("Running hard sphere hexane simulation");
        System.out.println(numMolecules + " molecules at density " + density);
        System.out.println("output data to " + filename);

        //spaces are now singletons; we can only have one instance, so we call
        // it with this method, not a "new" thing.
        TestHexaneHighMedium sim = new TestHexaneHighMedium(Space3D.getInstance(), numMolecules);

        PrimitiveHexane primitive = (PrimitiveHexane)sim.lattice.getPrimitive();
        
        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        waveVectorFactory = new WaveVectorFactorySimple(primitive);
        
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setPhase(sim.phase);

        IntervalActionAdapter fooAdapter = new IntervalActionAdapter(
                meterNormalMode);
        fooAdapter.setActionInterval(numMolecules);
        sim.integrator.addListener(fooAdapter);
        
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        meterNormalMode.reset();
        sim.getController().reset();
        
        sim.activityIntegrate.setMaxSteps(nSteps);      
        sim.getController().actionPerformed();
 
        System.out.println("Calculating all the fun eigenstuff");
        
        WriteS sWriter = new WriteS();
        sWriter.setFilename(filename);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setOverwrite(true);
        sWriter.actionPerformed();
        
        System.out.println(filename + " is finished running");
        System.out.println(nSteps + "  steps");
    }

}


