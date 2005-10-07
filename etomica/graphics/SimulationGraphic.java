package etomica.graphics;
import java.util.Iterator;
import java.util.LinkedList;

import javax.swing.JFrame;
import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.AtomPositionDefinition;
import etomica.integrator.Integrator;
import etomica.integrator.IntervalActionAdapter;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Polyhedron;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.simulation.SimulationContainer;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * General class for graphical presentation of the elements of a molecular simulation.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/26/02 (DAK) modified makeAndDisplayFrame method to return the frame
  * 09/13/02 (DAK) added blockDefaultLayout method.
  * 10/21/02 (DAK) added static method to set EtomicaTheme
  * 09/02/03 (DAK) setting Default.DO_SLEEP in constructor
  */
public class SimulationGraphic implements SimulationContainer, java.io.Serializable {
    
    static {
        try {
            javax.swing.plaf.metal.MetalLookAndFeel.setCurrentTheme(new EtomicaTheme());
//            javax.swing.plaf.metal.MetalLookAndFeel.setCurrentTheme(new BlueRoseTheme());
            javax.swing.UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
//            javax.swing.UIManager.setLookAndFeel(javax.swing.UIManager.getSystemLookAndFeelClassName());
        } catch(Exception e) {}
    }
    
    private SimulationPanel simulationPanel;
    
    public SimulationGraphic(Simulation simulation) {
        this.simulation = simulation;
        DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(simulation);
        add(controlPanel);
        setupDisplayPhase();
    }
    
    public Simulation getSimulation() {return simulation;}
    
    public final LinkedList displayList() { return displayList;}
    public final LinkedList deviceList() { return deviceList; }
        
    public void repaint() {
        
        
    }
    /**
     * A visual display of the simulation via a JPanel.
     */
     public SimulationPanel panel() {
        if(simulationPanel == null) simulationPanel = new SimulationPanel();
        return simulationPanel;
     }
     
     private void setupDisplayPhase() {
         Controller controller = simulation.getController();
         Action[] activities = controller.getPendingActions();
         LinkedList phaseList = new LinkedList();
         LinkedList integratorList = new LinkedList(); 
         for (int i=0; i<activities.length; i++) {
             if (activities[i] instanceof ActivityIntegrate) {
                 Integrator integrator = ((ActivityIntegrate)activities[i]).getIntegrator();
                 if (integratorList.contains(integrator)) continue;
                 integratorList.add(integrator);
                 Phase[] phases = integrator.getPhase();
                 for (int j=0; j<phases.length; j++) {
                     if (phaseList.contains(phases[j])) continue;
                     phaseList.add(phases[j]);
                     DisplayPhase display = new DisplayPhase(phases[j]);
                     add(display);
                     // check if 2D?
                     if (phases[j].space().D() != 3) {
                         integrator.addListener(new IntervalActionAdapter(display));
                     }
                 }
             }             
         }
     }

     public void add(Display display) {
         final java.awt.Component component = display.graphic(null);
         if(component == null) return; //display is not graphic
         if(display instanceof DisplayBox || display instanceof DisplayBoxesCAE) {
             final java.awt.GridBagConstraints gbcBox = new java.awt.GridBagConstraints();
             gbcBox.gridx = 0;
             panel().displayBoxPanel.add(component, gbcBox);
         }
         else {
             panel().displayPanel.add(display.getLabel(),component);
             //add a listener to update the tab label if the name of the display changes
             display.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
                 public void propertyChange(java.beans.PropertyChangeEvent evt) {
                     if(evt.getPropertyName().equals("label")) {
                         int idx = panel().displayPanel.indexOfComponent(component);
                         panel().displayPanel.setTitleAt(idx,evt.getNewValue().toString());
                     }
                 }
             });
         }
         displayList.add(display);
     }

     public void remove(Display display) {
         final java.awt.Component component = display.graphic(null);
         if(component == null) return; //display is not graphic
         if(display instanceof DisplayBox || display instanceof DisplayBoxesCAE) {
             panel().displayBoxPanel.remove(component);
         }
         else {
             panel().displayPanel.remove(component);
         }
         displayList.remove(display);
     }
     /**
      * Adds displays graphic to the simulation display pane
      */
     public void add(Device device) {
         java.awt.Component component = device.graphic(null);
         if(device instanceof DeviceTable) {
             panel().displayPanel.add(component);
         }
         else {
             final java.awt.GridBagConstraints gbc = new java.awt.GridBagConstraints();
             gbc.gridx = 0;
             panel().devicePanel.add(component,gbc);
         }
         deviceList.add(device);
     }
     
    public final JFrame makeAndDisplayFrame() {
        return makeAndDisplayFrame(panel());
    }
    
    public static JFrame makeAndDisplayFrame(JPanel panel) {
        JFrame f = new JFrame();
        f.setSize(700,500);
        f.getContentPane().add(panel);
        f.pack();
        f.show();
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
        return f;
    }
    
    public static final java.awt.event.WindowAdapter WINDOW_CLOSER 
        = new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        };
        
    public DisplayPhase getDisplayPhase(Phase phase) {
        Iterator iterator = displayList.iterator();
        while(iterator.hasNext()) {
            Object display = iterator.next();
            if(display instanceof DisplayPhase) {
                if(((DisplayPhase)display).getPhase() == phase) return (DisplayPhase)display;
            }
        }
        return null;
    }
    private final Simulation simulation;
    private final LinkedList displayList = new LinkedList();
    private final LinkedList deviceList = new LinkedList();

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
//        etomica.simulation.prototypes.SwMd2D sim = new etomica.simulation.prototypes.SwMd2D();
//        etomica.simulation.prototypes.LjMd2D sim = new etomica.simulation.prototypes.LjMd2D();
        etomica.simulation.prototypes.HsMc2d sim = new etomica.simulation.prototypes.HsMc2d();
//          etomica.simulation.prototypes.SWMD3D sim = new etomica.simulation.prototypes.SWMD3D();
//      etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D();
//      final etomica.simulation.prototypes.HSMD3DNoNbr sim = new etomica.simulation.prototypes.HSMD3DNoNbr();
//      etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
//        etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
//        etomica.simulation.prototypes.HSMD2D_atomNbr sim = new etomica.simulation.prototypes.HSMD2D_atomNbr();
//        etomica.simulation.prototypes.HSMD2D_noNbr sim = new etomica.simulation.prototypes.HSMD2D_noNbr();
//        etomica.simulation.prototypes.GEMCWithRotation sim = new etomica.simulation.prototypes.GEMCWithRotation();
      sim.getDefaults().doSleep = false;
      sim.getDefaults().ignoreOverlap = true;
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DeviceNSelector nSelector = new DeviceNSelector(sim,sim.phase.getAgent(sim.species));
        simGraphic.add(nSelector);
        
//        AtomFilterInPolytope filter = new AtomFilterInPolytope(sim.phase.boundary().getShape());
//        MyFilter filter = new MyFilter((Polyhedron)sim.phase.boundary().getShape());
//        PhaseDeleteMolecules deleter = new PhaseDeleteMolecules(filter);
        //positionDefinition shifts atom to same origin as polytope
//        AtomPositionDefinition position = new AtomPositionDefinition() {
//            public Vector position(Atom a) {
//                Vector3D r = (Vector3D)sim.phase.boundary().dimensions().clone();
//                r.TE(-0.5);
//                r.PE(a.coord.position());
//                return r;
//            }
//        };
//        deleter.setPhase(sim.phase);
//        filter.setPositionDefinition(position);
//        DeviceButton deleteButton = new DeviceButton(sim.getController(),deleter);
//        simGraphic.add(deleteButton);
        simGraphic.makeAndDisplayFrame();
        ColorSchemeByType.setColor(sim.species.getFactory().getType(), java.awt.Color.red);
//        ColorSchemeByType.setColor(sim.species2, java.awt.Color.blue);
        simGraphic.panel().setBackground(java.awt.Color.yellow);
        Plane plane = new Plane();
        plane.setThreePoints(new Vector3D(1,1,1), new Vector3D(2,2,2), new Vector3D(4,5,1));
        
    }//end of main
    
    private static class MyFilter implements AtomFilter, java.io.Serializable {
        Polyhedron polyhedron;
        AtomPositionDefinition positionDefinition;
        public MyFilter(Polyhedron shape) {
            this.polyhedron = shape;
        }
        public void setPositionDefinition(AtomPositionDefinition def) {
            positionDefinition = def;
        }
        public boolean accept(Atom atom) {
            Vector r = positionDefinition.position(atom);
            if(!polyhedron.contains(r)) return false;
            if(polyhedron.distanceTo(r) < 0.5) return false;
            return true;
            
        }
        
    }
}


