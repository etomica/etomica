package etomica.modules.crystalviewer;
import java.awt.GridLayout;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

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

    // minimum and maximum value for the slider.
    private int minimumPosition = -10;
    private int maximumPosition = 10;

    // Input widget unique IDs 
	private final int MILLER_INDEX_H = 0;
	private final int MILLER_INDEX_K = 1;
	private final int MILLER_INDEX_L = 2;
	private final int PLANE_SELECTION_BOX = 3;
	private final int POSITION_SLIDER = 4;

	private final int SLIDER_DECIMAL_PLACES = 2;

	private final double PLANE_TOLERANCE = 5.0e-3;
	

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

        ModifierLatticePlane modifier;
        latticePlane.setOrigin(new Vector3D());
        // Miller i indices
        boxA = new DeviceBox();
        modifier = new ModifierLatticePlane(MILLER_INDEX_H);
        modifier.setLabel("h");
        boxA.setModifier(modifier);
        boxA.setInteger(true);
        // Miller j indices
        boxB = new DeviceBox();
        modifier = new ModifierLatticePlane(MILLER_INDEX_K);
        modifier.setLabel("k");
        boxB.setModifier(modifier);
        boxB.setInteger(true);
        // Miller k indices
        boxC = new DeviceBox();
        modifier = new ModifierLatticePlane(MILLER_INDEX_L);
        modifier.setLabel("l");
        boxC.setModifier(modifier);
        boxC.setInteger(true);
        // position text box
        boxD = new DeviceBox();
        modifier = new ModifierLatticePlane(PLANE_SELECTION_BOX);
        modifier.setLabel("Plane Selection");
        boxD.setModifier(modifier);
        boxD.setInteger(false);
        boxD.setPrecision(2);
        // not sure
        boxH = new DeviceBox();
        boxH.setModifier(new ModifierLatticePlane(5));//4th index for hexagonal lattices
        boxH.setEditable(false);
        boxH.setInteger(true);

        millerPanel = new JPanel(new GridLayout(1,0));
        TitledBorder millerBorder = new TitledBorder("Miller Indices");
        millerBorder.setTitleJustification(TitledBorder.CENTER);
        millerPanel.setBorder(millerBorder);
        millerPanel.add(boxA.graphic());
        millerPanel.add(boxB.graphic());
        millerPanel.add(boxC.graphic());
        
        boxA.doUpdate();
        boxB.doUpdate();
        boxC.doUpdate();
        boxD.doUpdate();

        positionSlider = new DeviceSlider(null, new ModifierLatticePlane(POSITION_SLIDER));
        positionSlider.setPrecision(SLIDER_DECIMAL_PLACES);
        positionSlider.setMinimum(minimumPosition);
        positionSlider.setMaximum(maximumPosition);
        positionSlider.setNMajor(4);
        positionSlider.getSlider().setValue(0);
        positionSlider.setLabel("Distance From Origin");
        positionSlider.setShowBorder(true);
        positionSlider.setPostAction(new Action() {
            public void actionPerformed() {
                boxD.doUpdate();
                display.repaint();
            }
        });

        JPanel distancePanel = new JPanel();
        TitledBorder distanceBorder = new TitledBorder("Position");
        distanceBorder.setTitleJustification(TitledBorder.CENTER);
        distancePanel.setBorder(distanceBorder);
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

        latticePlane.setTolerance(PLANE_TOLERANCE);

        update();
    }

    public LatticePlane latticePlane() {return latticePlane;}

    public void update() {
//        if(latticePlane.isPrimitiveHexagonal()) {
//            millerPanel.add(boxH.graphic(), 2);
//            millerPanel.revalidate();
//        } else {
//            millerPanel.remove(boxH.graphic());
//            millerPanel.revalidate();
//        }
        positionSlider.doUpdate();
        display.repaint();
    }
    
    private class ModifierLatticePlane implements Modifier {
        
        private final int index;
        private String label = "modifier";
        
        public String getLabel() {
            return label;
        }
        public void setLabel(String newLabel) {
        	label = newLabel;
        }
        ModifierLatticePlane(int index) {
            this.index = index;
        }
        public void setValue(double t) {
            if(latticePlane == null) return;
            switch(index) {
                default:
                case MILLER_INDEX_H: 
                case MILLER_INDEX_K: 
                case MILLER_INDEX_L: 
                	latticePlane.setMillerIndex(index,(int)t);
                	break;
                case PLANE_SELECTION_BOX:
                	latticePlane.setPosition(t);
                	break;
                case POSITION_SLIDER:
                	latticePlane.setSpacePosition(t);
                	break;
                case 5: break; //4th hexagonal index
            }
            update();
        }
        public double getValue() {
            switch(index) {
                default:
                case MILLER_INDEX_H: 
                case MILLER_INDEX_K: 
                case MILLER_INDEX_L: return latticePlane.getMillerIndex(index);
                case PLANE_SELECTION_BOX: return latticePlane.getPosition();
                case POSITION_SLIDER: return latticePlane.getSpacePosition();
                case 5: return -(latticePlane.getMillerIndex(MILLER_INDEX_H) +
                		         latticePlane.getMillerIndex(MILLER_INDEX_K));
            }
        }
        public etomica.units.Dimension getDimension() {return Null.DIMENSION;}
    }
}