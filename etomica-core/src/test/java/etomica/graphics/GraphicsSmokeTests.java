package etomica.graphics;

import etomica.action.SimulationRestart;
import etomica.simulation.prototypes.HSMD3D;
import org.junit.After;
import org.junit.Test;

import javax.swing.*;

public class GraphicsSmokeTests {
    private JFrame frame;

    @After
    public void tearDown() {
        if(frame != null) {
            frame.dispose();
        }
    }

    @Test
    public void testHSMD3D() {
        HSMD3D sim = new HSMD3D();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "test", sim.getSpace(), sim.getController());
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
