package etomica.beans;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.*;
import java.beans.*;
import java.util.*;
import java.awt.*;
import java.io.*;
import etomica.*;
//import com.symantec.itools.vcafe.openapi.*;

public class PhaseCustomizer extends JTabbedPane implements Customizer
{
    etomica.Phase phase;
    BoundaryEditor boundaryEditor;
    Space.Boundary oldBoundary;
    JRadioButton[] boundaryButtons;
    DisplayPhase displayPhase;
    
    static Class[] meterClasses;
    
    public PhaseCustomizer()
    {
        //BOUNDARY PANEL
        boundaryEditor = new BoundaryEditor(this);
        String[] boundaries = boundaryEditor.boundaries();
        
        JPanel boundaryPane = new JPanel(new GridLayout(0,1));
        ButtonGroup boundaryButtonGroup = new ButtonGroup();
        boundaryButtons = new JRadioButton[boundaries.length];
        for(int i=0; i<boundaries.length; i++) {
            boundaryButtons[i] = new JRadioButton(boundaries[i],false);
            boundaryButtonGroup.add(boundaryButtons[i]);
            boundaryPane.add(boundaryButtons[i]);
            boundaryButtons[i].addActionListener(new BoundaryButtonListener());
        }
        
        //CONFIGURATION PANEL
        JPanel configurationPane = new JPanel();
        displayPhase = new DisplayPhase();
        configurationPane.add(displayPhase.graphic(null));
        displayPhase.addDisplayPhaseListener(new MoleculeMover());
        
        //METER PANEL
        JPanel meterPane = new JPanel(new GridLayout(0,2));
        for(int i=0; i<meterClasses.length; i++) {
            String name = meterClasses[i].getName();
            int idx = name.indexOf("r");  //strip off etomica.Meter prefix
            name = name.substring(idx+1);
            meterPane.add(new JCheckBox(name));
        }
        
           //workaround for JTabbedPane bug in JDK 1.2
        addChangeListener(
           new ChangeListener() {
               public void stateChanged(ChangeEvent event) {
                  validate();
               }
           });
        addTab("Configuration", configurationPane);
        addTab("Meters", meterPane);
        addTab("Boundary",boundaryPane);
        
        validate();
	}
	
	public void setObject(Object obj) {
	    phase = (Phase)obj;
	    
	    oldBoundary = phase.getBoundary();
	    for(int i=0; i<boundaryButtons.length; i++) {
	        if(boundaryButtons[i].getText().equals(oldBoundary.type().toString())) {
	            boundaryButtons[i].setSelected(true);
	            break;
	        }
	    }
	    displayPhase.setPhase(phase);
	    validate();
	    repaint();
	}
	
	public java.awt.Dimension getPreferredSize() {
	    return new java.awt.Dimension(400,400);
	}
	
	private class BoundaryButtonListener implements ActionListener, java.io.Serializable {
	    
	    public void actionPerformed(ActionEvent evt) {
	        //Action for selection of a new boundary via radio button
	        JToggleButton button = (JToggleButton)evt.getSource();
	        boundaryEditor.setAsText(button.getText());
	        phase.setBoundary((Space.Boundary)boundaryEditor.getValue());
	        PhaseCustomizer.this.firePropertyChange(
	                    "boundary", oldBoundary, phase.getBoundary());
	        oldBoundary = phase.getBoundary();
	    }//end of actionPerformed
	}//end of BoundaryButtonListener
	
	/**
	 * Override superclass method to make public so inner classes can fire events
	 */
	public void firePropertyChange(String property, Object oldValue, Object newValue) {
	    super.firePropertyChange(property, oldValue, newValue);
	}
	
	/**
	 * Listener to handle molecule repositioning in Customizer's configuration frame
	 */
	private class MoleculeMover implements DisplayPhaseListener, java.io.Serializable {
	    public void displayPhaseAction(DisplayPhaseEvent dpe) {
	        Molecule molecule = dpe.getMolecule();
	        if(molecule == null) return;
	        molecule.translateTo(dpe.getPoint());
	        ((Display)dpe.getSource()).repaint();
	        PhaseCustomizer.this.firePropertyChange(
	                    "moleculePosition", null, phase.getMoleculePosition());//phase.firstAtom().r);
	    }
	}
	/**
	 * Examines default class directory and locates all classes with names beginning with "Meter",
	 * then creates an array of Class files for each of these meter classes, excluding those that
	 * are not meant to be instantiated (e.g., MeterAbstract).
	 */
    static {
	    File dir = new File(etomica.Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
	            return name.startsWith("Meter") 
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.equals("Meter.class") 
	                && !name.equals("MeterFunction.class")
	                && !name.equals("MeterAbstract.class")
	                && !name.equals("MeterMultiFunction.class")
	                && !name.equals("MeterLocalDensity.class")
	                && !name.equals("MeterRatio.class");}
	        });
//	    java.util.Arrays.sort(files);  //won't autojar with this included
	    meterClasses = new Class[files.length];
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        meterClasses[i] = null;
	        try{
	           meterClasses[i] = Class.forName("etomica."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }
	}
	    
	
	public static void main(String[] args) {
	    File dir = new File(Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
	            return name.startsWith("Meter") && name.endsWith("class");}
	        });
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        Class cls = null;
	        try{
	           cls = Class.forName("etomica."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	        System.out.println(files[i]+"  "+cls.getName());
	    }
	}
}