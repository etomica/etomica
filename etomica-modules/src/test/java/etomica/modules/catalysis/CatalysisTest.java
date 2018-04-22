package etomica.modules.catalysis;

import etomica.action.ActionIntegrate;
import etomica.data.DataPumpListener;
import etomica.space3d.Space3D;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CatalysisTest {
    Catalysis sim;
    MeterDensityCO densityCO;
    MeterDensityO2 densityO2;
    MeterDensityCO2 densityCO2;
    ActionIntegrate actionIntegrate;

    @BeforeEach
    void setUp() {
        sim = new Catalysis(Space3D.getInstance(), 20);
        densityCO = new MeterDensityCO(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        densityO2 = new MeterDensityO2(sim.box, sim.speciesO, sim.interactionTracker.getAgentManager());
        densityCO2 = new MeterDensityCO2(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        actionIntegrate = new ActionIntegrate(sim.integrator);
        actionIntegrate.setMaxSteps(3000);
    }

    @Test
    void testDensities() {
        actionIntegrate.actionPerformed();
        assertEquals(4.1e-5, densityCO.getDataAsScalar(), 1e-5);
    }
}