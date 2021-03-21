package etomica.mu;

import etomica.graphics.SimulationGraphic;
import etomica.modules.mu.MuFasterer;
import etomica.modules.mu.MuGraphicFasterer;
import etomica.space.Space;
import org.assertj.swing.core.BasicRobot;
import org.assertj.swing.edt.GuiActionRunner;
import org.assertj.swing.fixture.FrameFixture;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import javax.swing.*;

public class MuGraphicTest {
    private FrameFixture frame;
    private MuFasterer sim;

    @BeforeEach
    public void setUp() {
        Space space = Space.getInstance(2);
        sim = new MuFasterer(space);
        JFrame f = GuiActionRunner.execute(() -> {
            MuGraphicFasterer swmdGraphic = new MuGraphicFasterer(sim);
            return SimulationGraphic.makeAndDisplayFrame
                    (swmdGraphic.getPanel(), "Chemical Potential");
        });
        frame = new FrameFixture(BasicRobot.robotWithCurrentAwtHierarchy(), f);
    }

    @Test
    public void test() throws InterruptedException {
        Assertions.assertDoesNotThrow(() -> {
            frame.button("runToggle").click();
            Thread.sleep(5000);
            frame.button("runToggle").click().requireText("Continue");
        });
        Assertions.assertTrue(sim.integrator.getStepCount() > 2);
    }
}
