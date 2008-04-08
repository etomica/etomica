package etomica.apps.junit.normalmode;

import junit.framework.TestCase;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.SimCalcSLJ;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
import etomica.space.Space;

/**
 * Tests calculating S for a simple (low T, low rho) LJ system.
 */
public class SimCalcSLJTest extends TestCase {

    /**
     * @param args
     */
    public static void main(String[] args) {
        new SimCalcSLJTest().testCalcS();
    }
    
    public void testCalcS() {

        // defaults
        int D = 3;
        int nA = 32;
        double density = 1.0;
        double temperature = .186;
        long simSteps = 50000;

        // construct simulation
        SimCalcSLJ sim = new SimCalcSLJ(Space.getInstance(D), nA, density, temperature);

        // set up initial configuration and save nominal positions
        Primitive primitive = sim.primitive;

        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        waveVectorFactory = new WaveVectorFactorySimple(primitive, sim.getSpace());
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);

        sim.integrator.addIntervalAction(meterNormalMode);
        sim.integrator.setActionInterval(meterNormalMode, nA);
        
        ((MCMoveStepTracker)sim.integrator.getMoveManager().getMCMoves()[0].getTracker()).setNoisyAdjustment(false);

        sim.activityIntegrate.setMaxSteps(simSteps/10);
        sim.getController().actionPerformed();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().reset();
        meterNormalMode.reset();

        WriteS sWriter = new WriteS(sim.getSpace());
        sWriter.setFilename("test");
        sWriter.setOverwrite(true);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setTemperature(temperature);
        
        sim.activityIntegrate.setMaxSteps(simSteps);
        sim.getController().actionPerformed();

        sWriter.actionPerformed();
        double A = sWriter.getLastA();
        
        assertTrue("harmonic reference free energy within expected limits "+A, Math.abs(A-45.3) < 0.3);
    }
}