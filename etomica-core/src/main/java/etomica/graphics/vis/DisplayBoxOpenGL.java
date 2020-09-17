package etomica.graphics.vis;

import com.jogamp.opengl.GLCapabilities;
import com.jogamp.opengl.GLProfile;
import com.jogamp.opengl.awt.GLCanvas;
import com.jogamp.opengl.util.Animator;
import com.jogamp.opengl.util.GLBuffers;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.Display;
import etomica.graphics.SimulationGraphic;
import etomica.math.geometry.LineSegment;
import etomica.simulation.prototypes.LJMD3D;
import etomica.simulation.prototypes.LJMD3DNbr;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import net.miginfocom.swing.MigLayout;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.nio.FloatBuffer;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;


public class DisplayBoxOpenGL extends Display {

    private final GLCanvas canvas;
    private final JPanel panel = new JPanel(new MigLayout());
    private final BoxRenderer renderer;
    private final Box box;
    private final Controller controller;
    private final Animator animator;

    public DisplayBoxOpenGL(Box box, Controller controller) {
        this.box = box;
        this.controller = controller;

        GLProfile profile = GLProfile.get(GLProfile.GL3);
        GLCapabilities caps = new GLCapabilities(profile);
        caps.setDepthBits(24);
        this.canvas = new GLCanvas(caps);
        this.panel.add(canvas);
        int boxEdgeCount = box.getBoundary().getShape().getEdges().length;
        this.renderer = new BoxRenderer(box.getLeafList().size(), boxEdgeCount);

        JCheckBox checkBox = new JCheckBox("Enable SSAO", true);
        checkBox.addItemListener(e -> {
            getRenderer().setEnableSSAO(e.getStateChange() != ItemEvent.DESELECTED);
        });
        this.panel.add(checkBox);

        ScheduledExecutorService scheduler = Executors.newScheduledThreadPool(1);
        ScheduledFuture<?> future = scheduler.scheduleAtFixedRate(() -> {
            controller.submitActionInterrupt(() -> {
                synchronized (renderer.positions) {
                    FloatBuffer pos = renderer.positions;
                    for (int i = 0; i < box.getLeafList().size(); i++) {
                        Vector vec = box.getLeafList().get(i).getPosition();
                        pos.put(i * 3, (float) vec.getX(0));
                        pos.put(i * 3 + 1, (float) vec.getX(1));
                        pos.put(i * 3 + 2, (float) vec.getX(2));
                    }
                }

                synchronized (renderer.boxVertices) {
                    FloatBuffer verts = renderer.boxVertices;
                    LineSegment[] lines = box.getBoundary().getShape().getEdges();
                    for (int i = 0; i < lines.length; i++) {
                        Vector v1 = lines[i].getVertices()[0];
                        Vector v2 = lines[i].getVertices()[1];

                        verts.put(i * 6 + 0, (float) v1.getX(0));
                        verts.put(i * 6 + 1, (float) v1.getX(1));
                        verts.put(i * 6 + 2, (float) v1.getX(2));
                        verts.put(i * 6 + 3, (float) v2.getX(0));
                        verts.put(i * 6 + 4, (float) v2.getX(1));
                        verts.put(i * 6 + 5, (float) v2.getX(2));
                    }
                }
            });
        }, 0, 33, TimeUnit.MILLISECONDS);
        scheduler.schedule(() -> future.cancel(false), 5, TimeUnit.SECONDS);

        canvas.addGLEventListener(renderer);
        canvas.addMouseMotionListener(renderer);
        canvas.addMouseWheelListener(renderer);
        canvas.setPreferredSize(new Dimension(800, 800));

        animator = new Animator(canvas);
        animator.setUpdateFPSFrames(60, System.err);

    }

    @Override
    public Component graphic() {
        return this.panel;
    }

    public void startAnimation() {
        animator.start();
    }

    public BoxRenderer getRenderer() {
        return this.renderer;
    }


    public static void main(String[] args) {
        LJMD3DNbr sim = new LJMD3DNbr();


        SwingUtilities.invokeLater(() -> {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim);
            graphic.remove(graphic.getDisplayBox(sim.box));
            DisplayBoxOpenGL display = new DisplayBoxOpenGL(sim.box, sim.getController());
            graphic.add(display);

            graphic.makeAndDisplayFrame();
            display.startAnimation();
        });
    }
}
