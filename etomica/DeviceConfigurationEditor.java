package etomica;
import javax.swing.*;
import java.awt.event.*;
import java.awt.Component;
import java.awt.GridLayout;

public class DeviceConfigurationEditor extends Device {
    
    public static boolean exists = false;
    private DisplayPhase display;
    private JInternalFrame buttonPane;
    String[] labels = new String[] {"No action", "Delete molecule", "Move molecule"};
    private final DisplayPhaseListener[] listeners = new DisplayPhaseListener[] {
        new NullListener(),
        new MoleculeDeleter(),
        new MoleculeMover()
    };
    DisplayPhaseListener currentListener;
    
    public DeviceConfigurationEditor(DisplayPhase display) {
        this(Simulation.instance, display);
    }
    public DeviceConfigurationEditor(Simulation sim, DisplayPhase display) {
        super(sim);
        this.display = display;
        exists = true;
        
        buttonPane = new JInternalFrame() {
            public void dispose() {
                parentSimulation().remove(this);
                parentSimulation().validate();
                parentSimulation().repaint();
                DeviceConfigurationEditor.exists = false;
                removeDisplayListeners();
                super.dispose();
            }//end of dispose
        };//end of inner subclass of JInternalFrame
        JPanel buttonPanel = new JPanel(new GridLayout(0,1));
        buttonPane.getContentPane().add(buttonPanel);
        buttonPane.setClosable(true);
        buttonPane.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        buttonPane.setTitle("Editor");
        ButtonGroup buttonGroup = new ButtonGroup();
        for(int i=0; i<labels.length; i++) {
            MyRadioButton button = new MyRadioButton(labels[i],i==0,listeners[i]);
            buttonGroup.add(button);
            buttonPanel.add(button);
            button.addActionListener(new ButtonListener());
        }
        
        currentListener = listeners[0];
        display.addDisplayPhaseListener(currentListener);
    }
    
    public void removeDisplayListeners() {
        display.removeDisplayPhaseListener(currentListener);
    }
    
    public Component graphic(Object obj) {return buttonPane;}
    
    /**
     * Simple extension of JRadioButton to permit it to be associated with an object
     */
    private class MyRadioButton extends JRadioButton {
        Object object;
        MyRadioButton(String label, boolean selected, Object obj) {
            super(label,selected);
            object = obj;
        }
    }
	private class ButtonListener implements ActionListener {
	    
	    public void actionPerformed(ActionEvent evt) {
	        MyRadioButton button = (MyRadioButton)evt.getSource();
	        display.removeDisplayPhaseListener(currentListener);
	        currentListener = (DisplayPhaseListener)button.object;
            display.addDisplayPhaseListener(currentListener);
	    }//end of actionPerformed
	}//end of ButtonListener
    
    //  Start of DisplayPhaseListeners //
    
	private class NullListener implements DisplayPhaseListener {
	    public void displayPhaseAction(DisplayPhaseEvent dpe) {}
	}
	
	private class MoleculeDeleter implements DisplayPhaseListener {
	    public void displayPhaseAction(DisplayPhaseEvent dpe) {
	        Atom molecule = dpe.getAtom();
	        if(molecule == null) return;
	        display.phase().removeMolecule(molecule);
	        display.repaint();
	        display.phase().integrator().initialize();
	    }
	}
    
	private class MoleculeMover implements DisplayPhaseListener {
	    public void displayPhaseAction(DisplayPhaseEvent dpe) {
	        Atom molecule = dpe.getAtom();
	        if(molecule == null) return;
	        molecule.coord.translateTo(dpe.getPoint());
	        display.repaint();
	        display.phase().integrator().initialize();
	    }
	}
    
    
}