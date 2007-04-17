package etomica.modules.crystalviewer;
import java.awt.GridLayout;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceCheckBox;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPhase;
import etomica.lattice.LatticePlane;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.space3d.Vector3D;
import etomica.units.Null;

/**
 * Class that produces a panel with controls that presence, orientation
 * and position of a plane that clips atoms from the display.  Permits
 * different surfaces of the lattice to be displayed.  Also controls
 * highlighting of atoms that are in the plane.
 */
 
 /* History
  * 09/07/02 (DAK) new
  */

public class ClipPlaneEditor {
    
    
    public ClipPlaneEditor(final LatticePlane latticePlane, final DisplayPhase display) {
    
        this.latticePlane = latticePlane;
        this.display = display;
        
        //box that toggles clipping display of atoms on one side of plane
        DeviceCheckBox clipToggle = new DeviceCheckBox("Clip", new ModifierBoolean() {
            public boolean getBoolean() {return display.getAtomFilter().equals(ClipPlaneEditor.this.latticePlane);}
            public void setBoolean(boolean b) {
                if(b) display.setAtomFilter(ClipPlaneEditor.this.latticePlane);
                else  display.setAtomFilter(null);
                display.repaint();
            }
        });
        
        //box that toggles highlighting of atoms in the plane
        colorScheme = new ColorSchemePlane(latticePlane,highlightColor,atomColor);
        display.setColorScheme(colorScheme);
        DeviceCheckBox highlightToggle = new DeviceCheckBox("Highlight", new ModifierBoolean() {
            public boolean getBoolean() {return colorScheme.getColorIn().equals(highlightColor);}
            public void setBoolean(boolean b) {
                if(b) colorScheme.setColorIn(highlightColor);
                else  colorScheme.setColorIn(atomColor);
                display.repaint();
            }
        });
        
        //box to toggle visibility of plane
        display.addDrawable(latticePlane);
        DeviceCheckBox showplaneToggle = new DeviceCheckBox("Show plane", new ModifierBoolean() {
            public boolean getBoolean() {return showPlane;}
            public void setBoolean(boolean b) {
                showPlane = b;
                if(b) display.addDrawable(latticePlane);
                else  display.removeDrawable(latticePlane);
            }
        });
        
        
        latticePlane.setOrigin(new Vector3D());
        boxA = new DeviceBox();
        boxA.setModifier(new ModifierLatticePlane(0));
        boxB = new DeviceBox();
        boxB.setModifier(new ModifierLatticePlane(1));
        boxC = new DeviceBox();
        boxC.setModifier(new ModifierLatticePlane(2));
        boxD = new DeviceBox();
        boxD.setModifier(new ModifierLatticePlane(3));
        boxH = new DeviceBox();
        boxH.setModifier(new ModifierLatticePlane(5));//4th index for hexagonal lattices
        boxH.setEditable(false);
        millerPanel = new JPanel(new GridLayout(1,0));
        millerPanel.setBorder(new javax.swing.border.TitledBorder("Miller indices"));
        millerPanel.add(boxA.graphic());
        millerPanel.add(boxB.graphic());
        millerPanel.add(boxC.graphic());
        
        boxA.setInteger(true);
        boxB.setInteger(true);
        boxC.setInteger(true);
        boxD.setInteger(true);
        boxH.setInteger(true);
        boxA.doUpdate();
        boxB.doUpdate();
        boxC.doUpdate();
        boxD.doUpdate();
        
        positionSlider = new DeviceSlider(null, new ModifierLatticePlane(4));
        positionSlider.setPrecision(1);
        positionSlider.setMinimum(-10);
        positionSlider.setMaximum(+10);
        positionSlider.setNMajor(4);
        positionSlider.getSlider().setValue(0);
        positionSlider.setPostAction(new Action() {
            public void actionPerformed() {
                boxD.doUpdate();
                display.repaint();
            }
        });
        
        JPanel distancePanel = new JPanel();
        distancePanel.setBorder(new javax.swing.border.TitledBorder("Position"));
        distancePanel.add(positionSlider.graphic());
        distancePanel.add(boxD.graphic());
        
        JPanel checkboxPanel = new JPanel(new java.awt.GridLayout(0,1));
        checkboxPanel.add(clipToggle.graphic());
        checkboxPanel.add(highlightToggle.graphic());
        checkboxPanel.add(showplaneToggle.graphic());
        
        int ix=0; 
        int iy=0;
        panel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc0 = new java.awt.GridBagConstraints();
        panel.add(checkboxPanel, gbc0);
        ix = 0;
        gbc0.gridx = ix; gbc0.gridy = ++iy;
        panel.add(millerPanel, gbc0);
        
        gbc0.gridx = ix; gbc0.gridy = ++iy;
        panel.add(distancePanel, gbc0);
        
        update();
    }
    public LatticePlane latticePlane() {return latticePlane;}
 
    public void update() {
        if(latticePlane.isPrimitiveHexagonal()) {
            millerPanel.add(boxH.graphic(), 2);
            millerPanel.revalidate();
        } else {
            millerPanel.remove(boxH.graphic());
            millerPanel.revalidate();
        }
        positionSlider.doUpdate();
    }
    
    public JPanel getPanel() {return panel;}
    
    private JPanel panel;
    protected DeviceBox boxA, boxB, boxC, boxD, boxH;
    protected LatticePlane latticePlane;
    private DeviceSlider positionSlider;
    
    protected DisplayPhase display;

    //highlight color for atoms when they are in the plane
    protected java.awt.Color highlightColor = java.awt.Color.yellow;
    protected java.awt.Color atomColor = java.awt.Color.red;
    protected ColorSchemePlane colorScheme;
    protected boolean showPlane = true;
    private JPanel millerPanel;

    
    private class ModifierLatticePlane implements Modifier {
        
        private final int index;
        public String getLabel() {
            return "modifier";
        }
        ModifierLatticePlane(int index) {
            this.index = index;
        }
        public void setValue(double t) {
            if(latticePlane == null) return;
            switch(index) {
                default:
                case 0: 
                case 1: 
                case 2: latticePlane.setMillerIndex(index,(int)t); break;
                case 3: latticePlane.setPosition((int)t); break;
                case 4: latticePlane.setSpacePosition(t); break;
                case 5: break; //4th hexagonal index
            }
            update();
        }
        public double getValue() {
            switch(index) {
                default:
                case 0: 
                case 1: 
                case 2: return latticePlane.getMillerIndex(index);
                case 3: return latticePlane.getPosition();
                case 4: return latticePlane.getSpacePosition();
                case 5: return -(latticePlane.getMillerIndex(0) + latticePlane.getMillerIndex(1));
            }
        }
        public etomica.units.Dimension getDimension() {return Null.DIMENSION;}
    }
}