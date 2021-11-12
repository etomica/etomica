package etomica.graphics;

import etomica.action.SimulationRestart;
import etomica.simulation.prototypes.HSMD3DFasterer;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;

import javax.swing.*;

public class GraphicsSmokeTests {
    private JFrame frame;

    @AfterEach
    public void tearDown() {
        if(frame != null) {
            frame.dispose();
        }
    }

    @Test
    public void testHSMD3D() {
        HSMD3DFasterer sim = new HSMD3DFasterer();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "test");
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);

        nSelector.setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(nSelector);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        frame = simGraphic.makeAndDisplayFrame();
    }
}
