package etomica.graphics;
import etomica.*;
import javax.swing.*;

/**
 * Class that permits selection and execution of the of any etomica class
 * that has a main method.
 * List of such classes is entered by hand in a static array (not generated automatically).
 *
 * @author David Kofke
 */
public class TestAll {
    
    static final String classes[] = {
        "ActionPointVolume", "AtomAction", "AtomIteratorSequential", "AtomPairIterator",
        "ColorSchemeBySpecies",
        "ColorSchemeHeteroAtoms", "ColorSchemeTemperature", "DeviceButton",
        "DeviceSlider", "DeviceTable", "DisplayBox", 
        "DisplayMeterHistogram", "DisplayTable", "DisplayToLog",
        "IntegratorGEMC", "IntegratorGear4", "MCMoveAtom", "MCMoveInsertDelete", "MCMovePointVolume",
        "MCMoveRotate", "MCMoveVolume", "MeterBondOrderParameterQ",
        "MeterPressureByVolumeChange", "MeterPressureHard", "MeterRDF",
        "Modulator", "PotentialAssociationCone",
        "PotentialFieldGravity", "PotentialRoughSphere", "Simulation",
        "SpeciesPistonCylinder", "SpeciesWalls"
    };
    static int idx = 0;
    static Class[] mainArgs = new Class[] {String[].class};
    
/*    public static void main(String[] args) {
        
        JFrame f = new JFrame();   //create a window
        JButton next = new JButton("Next");
        JButton finish = new JButton("Finish");
        JButton runSelected = new JButton("Run");
        JPanel panel = new JPanel(new java.awt.GridLayout(0,2));
        final ButtonGroup group = new ButtonGroup();
        for(int i=0; i<classes.length; i++) {
            JRadioButton radio = new javax.swing.JRadioButton(classes[i]);
            group.add(radio);
            panel.add(radio);
        }
        
        f.getContentPane().setLayout(new java.awt.FlowLayout());
//        f.getContentPane().add(next);
//        f.getContentPane().add(finish);
        f.getContentPane().add(runSelected);
        f.getContentPane().add(panel);
        
        next.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {doNext();}
        });
        finish.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {finish();}
        });
        runSelected.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                for(java.util.Enumeration e=group.getElements(); e.hasMoreElements(); ) {
                    JRadioButton button = (JRadioButton)e.nextElement();
                    if(button.isSelected()) {
                        doSim(button.getText());
                        break;
                    }
                }
                System.out.println();
            }
        });

        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
  */  
    public static void finish() {System.exit(0);}
    
    public static void doSim(String name) {
        String nextClass = "etomica."+name;
        System.out.println("Simulating "+nextClass);
        try {
            Class c = Class.forName(nextClass);
            java.lang.reflect.Method m = c.getMethod("main",mainArgs);
            m.invoke(null, new Object[] {new String[] {}});
        } 
        catch(ClassNotFoundException ex) {System.out.println("Class not available");}
        catch(NoSuchMethodException ex) {System.out.println("Method not available");}
        catch(java.lang.reflect.InvocationTargetException ex) {System.out.println("InvocationTargetException");}
        catch(IllegalAccessException ex) {System.out.println("Illegal access exception");}
        idx++;
    }     
    
    public static void doNext() {
        if(idx >= classes.length) finish();
        doSim(classes[idx]);
    }     
    
    
    public static Phase setupTestPhase(int nMolecules) {
        Simulation.instance = new Simulation();
        SpeciesSpheres species1 = new SpeciesSpheres(nMolecules);
        Phase phase = new Phase();
        Simulation.instance.elementCoordinator.go();
        return phase;
    }
}