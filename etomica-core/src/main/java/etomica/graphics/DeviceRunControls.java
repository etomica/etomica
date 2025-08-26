package etomica.graphics;

import etomica.action.controller.Controller;

import javax.swing.*;
import java.awt.*;

public class DeviceRunControls extends Device {
    private final JButton button;
    private boolean pauseDisplayed;

    public DeviceRunControls(Controller controller) {
        super(controller);
        this.button = new JButton(" Start ");
        this.button.setName("runToggle");
        this.pauseDisplayed = false;

        this.button.addActionListener(e -> {
            this.getController().start();
            if (controller.isPaused() == pauseDisplayed) {
                return;
            }
            this.controller.toggle().whenComplete((res, ex) -> {
                SwingUtilities.invokeLater(() -> {
                    if (this.controller.isPaused()) {
                        this.button.setText("Continue");
                        this.pauseDisplayed = false;
                    } else {
                        this.button.setText(" Pause ");
                        this.pauseDisplayed = true;
                    }
                });
            });
        });
    }

    public void reset() {
        SwingUtilities.invokeLater(() -> {
            this.button.setText(" Start ");
            this.pauseDisplayed = false;
        });
    }

    @Override
    public Component graphic() {
        return button;
    }
}
