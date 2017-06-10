/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.crystalviewer;
import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.simulation.Simulation;
import etomica.data.types.DataDouble;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceCheckBox;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
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
    protected DeviceBox boxH, boxK, boxL;
    protected LatticePlane latticePlane;
    private DeviceSlider positionSlider;
    private DisplayTextBox distanceDisplay;
    
    protected DisplayBox display;

    //highlight color for atoms when they are in the plane
    protected java.awt.Color highlightColor = java.awt.Color.yellow;
    protected java.awt.Color atomColor = java.awt.Color.red;
    protected ColorSchemePlane colorScheme;
    protected boolean showPlane = true;
    protected boolean showBoundary = true;
    private JPanel millerPanel;

    // minimum and maximum value for the slider.
    private int minimumPosition = -15;
    private int maximumPosition = 15;

    // Input widget unique IDs 
	private final int MILLER_INDEX_H = 0;
	private final int MILLER_INDEX_K = 1;
	private final int MILLER_INDEX_L = 2;
	private final int PLANE_SELECTION_SLIDER = 3;

	private final int SLIDER_DECIMAL_PLACES = 1;
	private final int DISTANCE_PRECISION    = 4;

	private final int MILLER_BOX_WIDTH  = 35;
	private final int MILLER_BOX_HEIGHT = 45;
	
	private final int MILLER_INDEX_MIN = -10;
	private final int MILLER_INDEX_MAX = 10;

    public ClipPlaneEditor(Simulation sim, final LatticePlane latticePlane, final DisplayBox display) {

        this.latticePlane = latticePlane;
        this.display = display;
        
        //box that toggles clipping display of atoms on one side of plane
        DeviceCheckBox clipToggle = new DeviceCheckBox("Clip", new ModifierBoolean() {
            public boolean getBoolean() {return display.getAtomFilter() != null;}
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
        ((DisplayBoxCanvasG3DSys)display.canvas).addPlane(latticePlane.getPlane());
        DeviceCheckBox showplaneToggle = new DeviceCheckBox("Show plane", new ModifierBoolean() {
            public boolean getBoolean() {return showPlane;}
            public void setBoolean(boolean b) {
                showPlane = b;
                if(b) ((DisplayBoxCanvasG3DSys)display.canvas).addPlane(latticePlane.getPlane());
                else  ((DisplayBoxCanvasG3DSys)display.canvas).removePlane(latticePlane.getPlane());
                display.canvas.repaint();
            }
        });

        // box to toggle crystal boundary display
        DeviceCheckBox crystalboundaryToggle = new DeviceCheckBox("Show boundary", new ModifierBoolean() {
            public boolean getBoolean() {return showBoundary;}
            public void setBoolean(boolean b) {
                showBoundary = b;
                display.setShowBoundary(b);
                display.repaint();
                }
        });

        ModifierLatticePlane modifier;
        latticePlane.setOrigin(new Vector3D());
        // Miller i indices
        boxH = new DeviceBox();
        modifier = new ModifierLatticePlane(MILLER_INDEX_H);
        modifier.setLabel("h");
        boxH.setModifier(modifier);
        boxH.setInteger(true);
        // Miller j indices
        boxK = new DeviceBox();
        modifier = new ModifierLatticePlane(MILLER_INDEX_K);
        modifier.setLabel("k");
        boxK.setModifier(modifier);
        boxK.setInteger(true);
        // Miller k indices
        boxL = new DeviceBox();
        modifier = new ModifierLatticePlane(MILLER_INDEX_L);
        modifier.setLabel("l");
        boxL.setModifier(modifier);
        boxL.setInteger(true);

        // Was not compiling in Java 1.4.  Component class (return
        // type of graphic() method) does not contain the setPreferredSize
        // method.  JComponent does.  I'm going to cast.
        JComponent j;
        j = (JComponent)boxH.graphic();
        j.setPreferredSize(new Dimension(MILLER_BOX_WIDTH, MILLER_BOX_HEIGHT));
        j = (JComponent)boxK.graphic();
        j.setPreferredSize(new Dimension(MILLER_BOX_WIDTH, MILLER_BOX_HEIGHT));
        j = (JComponent)boxL.graphic();
        j.setPreferredSize(new Dimension(MILLER_BOX_WIDTH, MILLER_BOX_HEIGHT));


        millerPanel = new JPanel(new GridLayout(1,0));
        TitledBorder millerBorder = new TitledBorder("Miller Indices");
        millerBorder.setTitleJustification(TitledBorder.CENTER);
        millerPanel.setBorder(millerBorder);
        millerPanel.add(boxH.graphic());
        millerPanel.add(boxK.graphic());
        millerPanel.add(boxL.graphic());

        boxH.doUpdate();
        boxK.doUpdate();
        boxL.doUpdate();

        distanceDisplay = new DisplayTextBox();
        distanceDisplay.setLabel("Distance from Origin");
        distanceDisplay.setPrecision(DISTANCE_PRECISION);

        positionSlider = new DeviceSlider(null, new ModifierLatticePlane(PLANE_SELECTION_SLIDER));
        positionSlider.setPrecision(SLIDER_DECIMAL_PLACES);
        // NEED TO SET MIN/MAX PLANE WHICH ARE DYNAMIC ...
        positionSlider.setMinimum(minimumPosition);
        positionSlider.setMaximum(maximumPosition);
        positionSlider.setNMajor(4);
        positionSlider.getSlider().setValue(0);
        positionSlider.setLabel("Plane Selection");
        positionSlider.setShowBorder(true);
        positionSlider.setShowValues(true);
        positionSlider.setEditValues(true);

        JPanel distancePanel = new JPanel();
        TitledBorder distanceBorder = new TitledBorder("Position");
        distanceBorder.setTitleJustification(TitledBorder.CENTER);
        distancePanel.setBorder(distanceBorder);
        distancePanel.add(positionSlider.graphic());
        distancePanel.add(distanceDisplay.graphic());

        JPanel checkboxPanel = new JPanel(new java.awt.GridLayout(2,2));
        TitledBorder checkboxBorder = new TitledBorder("Display Features");
        checkboxBorder.setTitleJustification(TitledBorder.CENTER);
        checkboxPanel.setBorder(checkboxBorder);
        checkboxPanel.add(clipToggle.graphic());
        checkboxPanel.add(highlightToggle.graphic());
        checkboxPanel.add(showplaneToggle.graphic());
        checkboxPanel.add(crystalboundaryToggle.graphic());

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

        // Add a componenet listner to the panel
        ClipPlaneEditorComponentListener listener =
        	                    new ClipPlaneEditorComponentListener();
        panel.addComponentListener(listener);

        update();
    }

    public LatticePlane latticePlane() {return latticePlane;}

    public void update() {
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
                	if ((int)t >= MILLER_INDEX_MIN && (int)t <= MILLER_INDEX_MAX) {
                		latticePlane.setMillerIndex(index,(int)t);
                	}
                	else if ((int)t < MILLER_INDEX_MIN) {
                	    latticePlane.setMillerIndex(index,MILLER_INDEX_MIN);
                	}
                	else {
                		latticePlane.setMillerIndex(index,MILLER_INDEX_MAX);
                	}
                	break;
                case PLANE_SELECTION_SLIDER:
                	latticePlane.setPosition(t);
                	break;
            }
        	DataDouble doubleD = new DataDouble();
        	doubleD.E(latticePlane.getSpacePosition());
            distanceDisplay.putData(doubleD);
            display.repaint();
        }

        public double getValue() {
            switch(index) {
                default:
                case MILLER_INDEX_H: 
                case MILLER_INDEX_K: 
                case MILLER_INDEX_L: return latticePlane.getMillerIndex(index);
                case PLANE_SELECTION_SLIDER: return latticePlane.getPosition();
            }
        }
        public etomica.units.Dimension getDimension() {return Null.DIMENSION;}
    }
    
    class ClipPlaneEditorComponentListener implements java.awt.event.ComponentListener {
        public void componentShown(java.awt.event.ComponentEvent ev) {
        	// When the screen is first shown, update the distance from
        	// origin.
        	DataDouble doubleD = new DataDouble();
        	doubleD.E(latticePlane.getSpacePosition());
            distanceDisplay.putData(doubleD);
        }
        public void componentMoved(java.awt.event.ComponentEvent ev) {
        }
        public void componentResized(java.awt.event.ComponentEvent ev) {
        }
        public void componentHidden(java.awt.event.ComponentEvent ev) {
        }    	
    }
}
