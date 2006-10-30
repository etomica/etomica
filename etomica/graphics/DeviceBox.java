//class includes a main method to demonstrate and test its use
package etomica.graphics;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.modifier.Modifier;
import etomica.modifier.ModifyAction;
import etomica.units.systems.UnitSystem;
import etomica.util.Constants;
import etomica.util.EnumeratedType;

/**
 * A simple device the permits editing of a single value via a textbox 
 * with an associated label.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/05/02 (DAK) new from DisplayBox
  * 09/22/02 (DAK) added integer flag to indicate display as decimal or integer
  * 09/24/02 (DAK) added key listener that increments/decrements with arrow keys, if integer
  * 09/27/02 (DAK) modified BoxListener so that it handles NumberFormatException appropriately
  *                added set/is Editable fields
  */
 
public class DeviceBox extends Device implements EtomicaElement, javax.swing.event.ChangeListener {
    
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel label;
    protected String labelString;
    private Constants.CompassDirection labelPosition = Constants.CompassDirection.NORTH;
    /**
     * Object for displaying the value as a text field
     */
    protected JTextField textField;
    /**
     * Displayed panel that holds the label and value
     * (not yet used; meant to implement to make lightweight display)
     */
    protected JPanel panel;
    /**
     * Modifier connecting the slider to the property
     */
    protected ModifyAction modifyAction;
    /**
     * Integer specifying the number of significant figures to be displayed.
     * Default is 4.
     */
    protected int precision;
    private LabelType labelType;
    private boolean integer = false;
    
    public DeviceBox() {
        super();
        panel = new JPanel(new java.awt.BorderLayout());
        label = null;
        labelString = null;
        textField = new JTextField("");
        textField.setEditable(true);
        panel.add(textField, java.awt.BorderLayout.CENTER);
        setLabelType(LabelType.STRING);
        unit = etomica.units.Null.UNIT;
        setPrecision(4);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple textbox editor of a single value");
        return info;
    }
    
    /** 
     * calls doUpdate method.  Implementation of ChangeListener interface.
     */
    public void stateChanged(javax.swing.event.ChangeEvent evt) {
        doUpdate();
    }
    
    /**
     * Updates the display of the box with the current value given by the modifier.
     */
    protected void doUpdate() {
        if(modifyAction == null) return;
        if(integer) {
            textField.setText(Integer.toString((int)unit.fromSim(modifyAction.getWrappedModifier().getValue())));
        }
        else {
            textField.setText(DisplayBox.format(unit.fromSim(modifyAction.getWrappedModifier().getValue()),precision));
        }
    }
    
    public void setEditable(boolean b) {textField.setEditable(b);}
    public boolean isEditable() {return textField.isEditable();}
    
    /**
     * Accessor method to set the physical units of the displayed value.
     * Text describing unit is used in label.
     */
    public void setUnit(etomica.units.Unit u) {
        unit = u;
        doUpdate();
        setLabel(labelString);
    }
    
    public java.awt.Component graphic(Object obj) {return panel;}
    
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {return precision;}

    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
        textField.setColumns(n+2);
        precision = n;
    }
    
    /**
     * Sets a flag indicating if the value should be displayed as an integer.
     */
    public void setInteger(boolean b) {integer = b;}
    public boolean isInteger() {return integer;}
    
    /**
     * Specifies the modifier that receives the edit.
     */
    public void setModifier(Modifier m) {
        if(m == null) throw new NullPointerException();
        if(unit == null) {
            setUnit(m.getDimension().getUnit(UnitSystem.SIM));
        }
        modifyAction = new ModifyAction(m);
        setLabel(labelString);
        doUpdate();

        BoxListener listener = new BoxListener();
        textField.addActionListener(listener);
        textField.addKeyListener(listener);
        textField.addFocusListener(listener);
    }
    
    /**
     * Accessor method for the modifier that receives the edit.
     */
    public Modifier getModifier() {
        return modifyAction == null ? null : modifyAction.getWrappedModifier();
    }
    
    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s) {
        labelString = s;
        if(labelType == LabelType.BORDER) {
            panel.setBorder(new javax.swing.border.TitledBorder(getLabel()));
        }
        else if(labelType == LabelType.STRING) {
            label.setText(getLabel());
            setLabelPosition(labelPosition);
        }
    }
    
    /**
     * @return the current value of the descriptive label.
     */
    public String getLabel() {
        if (labelString != null) {
            return labelString;
        }
        if (modifyAction == null) {
            return "";
        }
        String symbolText = unit.symbol();
        return modifyAction.getLabel()+
                ((symbolText.length() > 0) ? " ("+symbolText+")" : "");
    }
    
    /**
     * Sets the label type to "border" or "string".  With border, the label
     * text is part of the border.  With "string", the label is inside the 
     * border.
     */
    public void setLabelType(LabelType newLabelType) {
        if (newLabelType == labelType) {
            return;
        }
        labelType = newLabelType;
        if(labelType == LabelType.BORDER) {
            panel.remove(label);
            label = null;
            panel.setBorder(new javax.swing.border.TitledBorder(getLabel()));
        }
        else if(labelType == LabelType.STRING) {
            label = new JLabel(getLabel());
            panel.setBorder(new javax.swing.border.EmptyBorder(2,2,2,2));
            setLabelPosition(labelPosition);
        }
    }

    public LabelType getLabelType() {
        return labelType;
    }

    public void setLabelPosition(Constants.CompassDirection position) {
        labelPosition = position;
        if(labelType != LabelType.STRING) return;
        panel.remove(label);
        panel.add(label,position.toString());//toString() returns the corresponding BorderLayout constant
        panel.revalidate();
        panel.repaint();
    }
    
    public Constants.CompassDirection getLabelPosition() {
        return labelPosition;
    }
        
    private class BoxListener extends java.awt.event.KeyAdapter implements java.awt.event.ActionListener, FocusListener {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            updateValue();
        }
        public void keyReleased(java.awt.event.KeyEvent evt) {
            if(!isInteger()) return;
            int key = evt.getKeyCode();
            int step = 0;
            if(key == KeyEvent.VK_UP || key == KeyEvent.VK_RIGHT) step = 1;
            else if(key == KeyEvent.VK_DOWN || key == KeyEvent.VK_LEFT) step = -1;
            if(step == 0) return;
            int i = Integer.parseInt(textField.getText()) + step;
            textField.setText(Integer.toString(i));
            updateValue();
        }
       
        public void focusGained(FocusEvent evt) {
        }
        public void focusLost(FocusEvent evt) {
            // user tabbed or clicked out of the text field
            updateValue();
       }
       
       /**
        * Take the value form the textfield and apply it to the simulation via 
        * the modifyAction.  Update the textbox accordingly. 
        */
       protected void updateValue() {
           double oldX = modifyAction.getWrappedModifier().getValue();
           double newX;
           try { 
               newX = unit.toSim(Double.parseDouble(textField.getText()));
           } catch(NumberFormatException ex) {//if bad format, just restore original value
               doUpdate();
               return;
           }
           if (newX != oldX) {
               modifyAction.setValueForAction(newX);
               try {
                   doAction(modifyAction);
               }
               catch (RuntimeException e) {
                   // if the modifier throws, we'll revert to the old value in doUpdate
               }
           }
           doUpdate();
       }
    }
    
    /**
     * Typed constant used to indicate the type of label to be used with the display.
     */
	public static class LabelType extends EnumeratedType {
        public LabelType(String label) {super(label);}       
        public static final LabelType BORDER = new LabelType("Border");
        public static final LabelType STRING = new LabelType("String");

        public static final LabelType[] choices() { 
            return new LabelType[] {BORDER,STRING};
        }
    }
    
}