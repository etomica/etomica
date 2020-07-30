package etomica.graphics;

import etomica.action.activity.Controller2;

import javax.swing.*;
import java.awt.*;

public class DeviceRunControls extends Device {
    private final Controller2 controller;
    private final JButton button;
    private boolean pauseDisplayed;

    public DeviceRunControls(Controller2 controller) {
        this.controller = controller;
        this.button = new JButton(" Start ");
        this.pauseDisplayed = false;

        this.button.addActionListener(e -> {
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
