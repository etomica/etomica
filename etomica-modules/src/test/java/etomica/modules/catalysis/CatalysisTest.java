package etomica.modules.catalysis;

import etomica.action.ActionIntegrate;
import etomica.space3d.Space3D;
import etomica.util.random.RandomMersenneTwister;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CatalysisTest {
    private static final double delta = 1e-9;
    Catalysis sim;
    MeterDensityCO densityCO;
    MeterDensityO2 densityO2;
    MeterDensityCO2 densityCO2;
    ActionIntegrate actionIntegrate;

    @BeforeEach
    void setUp() {
        sim = new Catalysis(Space3D.getInstance(), 20, new RandomMersenneTwister(1));
        densityCO = new MeterDensityCO(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        densityO2 = new MeterDensityO2(sim.box, sim.speciesO, sim.interactionTracker.getAgentManager());
        densityCO2 = new MeterDensityCO2(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        actionIntegrate = new ActionIntegrate(sim.integrator);
        actionIntegrate.setMaxSteps(10000);
    }

    @Test
    void testDensities() {
        assertDoesNotThrow(() -> actionIntegrate.actionPerformed());
        assertAll(
                () -> assertEquals(3.5144282273534565e-5, densityCO.getDataAsScalar(), delta, "CO density"),
                () -> assertEquals(4.1001662652457E-5, densityO2.getDataAsScalar(), delta, "O2 density"),
                () -> assertEquals(1.1714760757844856E-6, densityCO2.getDataAsScalar(), delta, "CO2 density")
        );

    }
}
