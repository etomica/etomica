package etomica.graphics.vis;

import com.jogamp.opengl.GLCapabilities;
import com.jogamp.opengl.GLProfile;
import com.jogamp.opengl.awt.GLCanvas;
import com.jogamp.opengl.util.Animator;
import com.jogamp.opengl.util.GLBuffers;
import etomica.action.controller.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.Display;
import etomica.graphics.SimulationGraphic;
import etomica.simulation.prototypes.LJMD3D;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

import javax.swing.SwingUtilities;
import java.awt.Component;
import java.awt.Dimension;
import java.nio.FloatBuffer;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;


public class DisplayBoxOpenGL extends Display {

    private final GLCanvas canvas;
    private final BoxRenderer renderer;
    private final Box box;
    private final Controller controller;

    public DisplayBoxOpenGL(Box box, Controller controller) {
        this.box = box;
        this.controller = controller;

        GLProfile profile = GLProfile.get(GLProfile.GL3);
        GLCapabilities caps = new GLCapabilities(profile);
        this.canvas = new GLCanvas(caps);
        this.renderer = new BoxRenderer(box.getLeafList().size());

        ScheduledExecutorService scheduler = Executors.newScheduledThreadPool(1);
        scheduler.scheduleAtFixedRate(() -> {
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
            });
        }, 0, 33, TimeUnit.MILLISECONDS);

        canvas.addGLEventListener(renderer);
        canvas.setPreferredSize(new Dimension(800, 800));

        Animator animator = new Animator(canvas);
        animator.setUpdateFPSFrames(60, System.err);
        animator.start();

    }

    @Override
    public Component graphic() {
        return this.canvas;
    }


    public static void main(String[] args) {
        LJMD3D sim = new LJMD3D();


        SwingUtilities.invokeLater(() -> {
            SimulationGraphic graphic = new SimulationGraphic(sim);
            graphic.remove(graphic.getDisplayBox(sim.box));
            graphic.add(new DisplayBoxOpenGL(sim.box, sim.getController()));

            graphic.makeAndDisplayFrame();
        });
    }
}
