package etomica;
import etomica.units.*;
import javax.swing.*;
import java.awt.event.*;
import java.awt.Component;
import java.awt.GridLayout;

public class DeviceUnitEditor extends Device {
    
    private Dimensioned element;
    private Prefix prefix = Prefix.NULL;
    private BaseUnit baseUnit = BaseUnit.Null.UNIT;
    private Dimension dimension;
    private JInternalFrame buttonPane;
    private Unit oldUnit;
    
//    private JComboBox prefixes = new JComboBox(Prefix.ALL);
    private JComboBox prefixes = new JComboBox();
    
    public DeviceUnitEditor(final Dimensioned element) {
        this(Simulation.instance, element);
    }
    public DeviceUnitEditor(Simulation sim, final Dimensioned element) {
        super(sim);
        this.element = element;
        
        oldUnit = element.getUnit();
        dimension = element.dimension();
        prefix = oldUnit.getPrefix();
  //      Class[] classes = getClasses(dimension);
        Class[] classes = BaseUnit.all(dimension);
        
        for(int i=0; i<Prefix.ALL.length; i++) {
            prefixes.addItem(Prefix.ALL[i]);
        }
        prefixes.setSelectedItem(prefix);
        prefixes.setKeySelectionManager(new KeySelector());
        prefixes.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                prefix = (Prefix)prefixes.getSelectedItem(); 
                element.setUnit(new Unit(prefix,baseUnit));
                buttonPane.repaint();
            }});
            
            
        buttonPane = new JInternalFrame() {
            public void dispose() {
         //       element.setUnit(new Unit(prefix,baseUnit));
                parentSimulation().remove(this);
                parentSimulation().validate();
                parentSimulation().repaint();
                super.dispose();
            }//end of dispose
        };//end of inner subclass of JInternalFrame
        JPanel buttonPanel = new JPanel(new GridLayout(0,1));
        buttonPane.getContentPane().add(buttonPanel);
        buttonPane.setClosable(true);
        buttonPane.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        buttonPane.setTitle("Units Editor");
        ButtonGroup buttonGroup = new ButtonGroup();
        buttonPanel.add(prefixes);
        for(int i=0; i<classes.length; i++) {
            StringBuffer name = new StringBuffer(classes[i].getName());
            name.delete(0,14);  //drop "etomica.units."
            MyRadioButton button = new MyRadioButton(name.toString(),oldUnit.getClass().getName().equals(classes[i].getName()),classes[i]);
            buttonGroup.add(button);
            buttonPanel.add(button);
            button.addActionListener(new ButtonListener());
        }
    }
    
    public Component graphic(Object obj) {return buttonPane;}
    
    /**
     * Simple extension of JRadioButton to permit it to be associated with an object
     */
    private class MyRadioButton extends JRadioButton {
        Class cls;
        MyRadioButton(String label, boolean selected, Class cls) {
            super(label,selected);
            this.cls = cls;
        }
    }
	private class ButtonListener implements ActionListener {
	    
	    public void actionPerformed(ActionEvent evt) {
	        MyRadioButton button = (MyRadioButton)evt.getSource();
	        try {
	            baseUnit = (BaseUnit)button.cls.newInstance();
	        }
	        catch(InstantiationException e) {System.out.println(e.toString()); System.exit(1);}
	        catch(IllegalAccessException e) {System.out.println(e.toString()); System.exit(1);}
            element.setUnit(new Unit(prefix,baseUnit));
	    }//end of actionPerformed
	}//end of ButtonListener
		
	private class KeySelector implements JComboBox.KeySelectionManager {
	    public int selectionForKey(char aKey, ComboBoxModel aModel) {
	        return Prefix.intKeySelect(aKey);
	    }
	}
}