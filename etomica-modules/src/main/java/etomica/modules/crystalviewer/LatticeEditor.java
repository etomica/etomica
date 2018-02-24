/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.crystalviewer;

import etomica.action.IAction;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceSlider;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Primitive;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Degree;
import etomica.units.dimensions.*;
import etomica.util.Arrays;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.util.HashMap;

/**
 * Class that produces a panel with controls that permit editing of
 * the structure of a lattice.  Permits changing of type primitive and
 * parameters (sizes and angles) for lattice vectors.
 */
public class LatticeEditor {
    
    public static int DEFAULT_SIZE = 6;
    protected int size;
    protected BravaisLattice currentLattice;
    private JPanel panel;
    protected Box box;
    protected ISpecies species;
    private DeviceBox[] angleBoxes, sizeBoxes;
    public JPanel anglePanel, sizePanel;
    public JPanel boxPanel;
    protected CrystalViewer viewer;
    protected final HashMap latticeNameHash;
    protected LatticeEditorBoxPropertyArray pvBox = null; 
    private final Space space;
    
    private final String[] fieldTitles = {"A", "B", "C",
                                          "Alpha", "Beta", "Gamma" };
    private final String[] fieldPrefix = {"size", "size", "size",
                                          "angle", "angle", "angle"};
    private final Color[] colorList = { Color.RED,
    		                            Color.GREEN,
    		                            Color.BLUE,
    		                            Color.CYAN,
    		                            Color.YELLOW,
    		                            new Color(255, 0, 255),
    		                            Color.WHITE,
    		                            Color.WHITE,
    		                            Color.WHITE,
    		                            Color.WHITE
                                      };

    private static final int UNIT_CELL_START_INDEX = 0;
    private static final int UNIT_CELL_END_INDEX = 2;
    private static final int ANGLE_START_INDEX = 3;
    private static final int ANGLE_END_INDEX = 5;

    public LatticeEditor(CrystalViewer viewer, BravaisLattice[] lattices,
    		             String[] latticeNames, Space _space) {
    
        this.viewer = viewer;
        this.space = _space;
        box = viewer.box;
        species = viewer.species;
        latticeNameHash = new HashMap();
        for (int i=0; i<lattices.length; i++) {
            latticeNameHash.put(lattices[i], latticeNames[i]);
        }
        final javax.swing.JComboBox selector = new javax.swing.JComboBox(lattices);
        selector.setRenderer(new javax.swing.plaf.basic.BasicComboBoxRenderer() {
            public Component getListCellRendererComponent(JList list, Object value, 
                    int index, boolean isSelected, boolean cellHasFocus) {
                //this uses lattice.toString(), which is bogus
                super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                //override the string with the action lattice name from the hash 
                String latticeName = (String)latticeNameHash.get(value);
                setText((latticeName == null) ? "" : latticeName);
                return this;
            }
        });
        currentLattice = lattices[0];
        selector.setSelectedIndex(0);
        
        sizeBoxes = new DeviceBox[0];
        angleBoxes = new DeviceBox[0];
        
        selector.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                currentLattice = (BravaisLattice)selector.getSelectedItem();
                changeBox();
                initBoxes();
                update();
            }
        });

        anglePanel = new JPanel(new GridLayout(1,3));
        GridLayout gl = new GridLayout(1,3);
        gl.setHgap(2);
        sizePanel = new JPanel(gl);
        sizePanel.setBorder(new TitledBorder("Cell Dimensions"));
        DeviceSlider nSlider = new DeviceSlider(null, new NModifier());
        nSlider.setPrecision(0);
        nSlider.setMinimum(1);
        nSlider.setMaximum(10);
        nSlider.setNMajor(9);
        nSlider.getSlider().setValue(DEFAULT_SIZE);
        nSlider.setBorderAlignment(TitledBorder.CENTER);
        nSlider.setLabel("Unit Cells Per Side");
        nSlider.setShowBorder(true);

        boxPanel = new JPanel(new GridLayout(2,1));
        boxPanel.add(sizePanel);
        TitledBorder bp = new TitledBorder("Primitive Vectors");
        bp.setTitleJustification(TitledBorder.CENTER);
        boxPanel.setBorder(bp);
       
        createBoxes();
        initBoxes();
        
        int ix = 0; int iy = 0;
        panel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc0 = new java.awt.GridBagConstraints();
        panel.add(selector, gbc0);
        gbc0.gridx = ix; gbc0.gridy = ++iy;
        panel.add(boxPanel, gbc0);
        gbc0.gridx = ix; gbc0.gridy = ++iy;
        panel.add(nSlider.graphic(), gbc0);
        
        update();

    }
    
    /**
     * Returns the currently selected crystal.
     */
    public BravaisLattice getCurrentLattice() {
        return currentLattice;
     }
    
    private void createBoxes() {

        angleBoxes = new DeviceBox[0];
        sizeBoxes = new DeviceBox[0];

        // Create a lattice editor property box array
        pvBox = new LatticeEditorBoxPropertyArray(this.fieldTitles, this.fieldPrefix);

    	// Unit cell size entry fields
        for(int box = UNIT_CELL_START_INDEX; box <= UNIT_CELL_END_INDEX; box++) {

            DeviceBox newBox = new DeviceBox();

            sizePanel.add(newBox.graphic());
            sizeBoxes = (DeviceBox[])Arrays.addObject(sizeBoxes, newBox);

            newBox.setPostAction(new IAction() {
                public void actionPerformed() {
                    update();
                }
            });
        }

        // Angle entry fields
    	for(int angle = ANGLE_START_INDEX; angle <= ANGLE_END_INDEX; angle++) {

    		DeviceBox newBox = new DeviceBox();

            newBox.setUnit(Degree.UNIT);
            anglePanel.add(newBox.graphic());
            anglePanel.setBorder(new TitledBorder("Angles"));
            angleBoxes = (DeviceBox[])Arrays.addObject(angleBoxes, newBox);

            newBox.setPostAction(new IAction() {
                public void actionPerformed() {
                    update();
                }
            });
    	}

        boxPanel.setLayout(new GridLayout(2,1));
        boxPanel.add(anglePanel);

    }

    private void initBoxes() {

    	// Unregister all value changed events
        for(int box = UNIT_CELL_START_INDEX; box <= UNIT_CELL_END_INDEX; box++) {
        	sizeBoxes[box].removeAllValueChangedListeners();
        }

        // Need to initialze box array here
        pvBox.initialize();

        ModifierGeneral modifier;

        // Unit cell size entry fields
        for(int box = UNIT_CELL_START_INDEX; box <= UNIT_CELL_END_INDEX; box++) {

            modifier = pvBox.getModifier(this.fieldTitles[box]);

            if (pvBox.hasWrite(this.fieldTitles[box]) == false) {
            	sizeBoxes[box-UNIT_CELL_START_INDEX].setLabel(fieldTitles[box]);
            	sizeBoxes[box-UNIT_CELL_START_INDEX].setEditable(false);
                // If there is no modifier for the box, the DeviceBox will
                // throw a NullPointerException when its modifier is set.
                // Discard the exception.
                try {
                	sizeBoxes[box-UNIT_CELL_START_INDEX].setModifier(modifier);
                }
                catch (NullPointerException ex) {}

                sizeBoxes[box-UNIT_CELL_START_INDEX].doUpdate();
            }
            else {
            	sizeBoxes[box-UNIT_CELL_START_INDEX].setEditable(true);
            	sizeBoxes[box-UNIT_CELL_START_INDEX].setModifier(modifier);
            	sizeBoxes[box-UNIT_CELL_START_INDEX].setBorderBackground
            	                (colorList[pvBox.getRootIndex(this.fieldTitles[box])]);
            }
        }

        // Angle entry fields
        for(int abox = ANGLE_START_INDEX; abox <= ANGLE_END_INDEX; abox++) {
            modifier = pvBox.getModifier(this.fieldTitles[abox]);

            if (pvBox.hasWrite(this.fieldTitles[abox]) == false) {
                angleBoxes[abox-ANGLE_START_INDEX].setLabel(fieldTitles[abox]);
                angleBoxes[abox-ANGLE_START_INDEX].setEditable(false);
                // If there is no modifier for the box, the DeviceBox will
                // throw a NullPointerException when its modifier is set.
                // Discard the exception.
                try {
                    angleBoxes[abox-ANGLE_START_INDEX].setModifier(modifier);
                }
                catch (NullPointerException ex) {}

                angleBoxes[abox-ANGLE_START_INDEX].doUpdate();
            }
            else {
        	    angleBoxes[abox-ANGLE_START_INDEX].setEditable(true);
                angleBoxes[abox-ANGLE_START_INDEX].setModifier(modifier);
            }
        }

        registerValueChangedEvents(sizeBoxes);
    }

    private void registerValueChangedEvents(DeviceBox[] boxes) {

    	int numBoxes = boxes.length;

    	for(int i = 0; i < numBoxes; i++) {
    		for(int j = 0; j <numBoxes; j++) {
    			if (i != j) {
    				boxes[i].addValueChangedListener(boxes[j]);
    			}
    		}
    	}
    }

    protected void update() {
        box.getBoundary().setBoxSize(Vector.of(size, size, size));

        int numAtoms = size*size*size;
        if (currentLattice instanceof BravaisLatticeCrystal) {
            numAtoms *= ((BravaisLatticeCrystal)currentLattice).getBasis().getScaledCoordinates().length;
        }
        box.setNMolecules(species, numAtoms);
        ConfigurationLattice config = new ConfigurationLattice(currentLattice, space);
        config.initializeCoordinates(box);
        viewer.update(currentLattice);
    }

    protected void changeBox() {

        Box oldBox = box;
        double[] boxSize = new double[]{10.0, 10.0, 10.0};
        if (oldBox != null) {
            viewer.sim.removeBox(oldBox);
        }

        box = viewer.sim.makeBox(new BoundaryDeformableLattice(currentLattice.getPrimitive(), boxSize));
        viewer.displayBox.setBox(box);
    }

    public JPanel getPanel() {return panel;}

    protected class NModifier implements Modifier {
        public void setValue(double d) {
            if ((int)d == size) return;
            if(d > 10) d = DEFAULT_SIZE;
            size = (int)d;
            update();
        }
        public double getValue() {
            return size;
        }
        public etomica.units.dimensions.Dimension getDimension() {return Quantity.DIMENSION;}
        public String getLabel() {
            return "a label";
        }
    } // end protected class NModifier

    /**
     * LatticeEditorBoxPropertyArray
     * @author rrassler
     *
     */
    private class LatticeEditorBoxPropertyArray {
    	private int numItems = 0;
    	private LatticeEditorBoxProperty[] boxes = null;
    	private int[] rootIndex = null;

    	public LatticeEditorBoxPropertyArray(String[] titles, String[] prefix) {

    		if (titles.length != prefix.length) {
    			this.numItems = 0;
    			return;
    		}
    		else {
    			this.numItems = titles.length;
    			rootIndex = new int[this.numItems];
    			rootIndex[0] = 0;
    		}

    		boxes = new LatticeEditorBoxProperty[numItems];

    		// Create property boxes
    		for(int i = 0; i < numItems; i++) {
			    boxes[i] = new LatticeEditorBoxProperty(titles[i], prefix[i], null);
    		}

    	}

    	public void initialize() {

            Primitive primitive = currentLattice.getPrimitive();
            Class c = primitive.getClass();
            //bits of this code are taken from Thinking in Java (1st edition), pages 708-713
            BeanInfo bi = null;
            try {
                bi = Introspector.getBeanInfo(c, java.lang.Object.class);
            }
            catch(IntrospectionException ex) {
                System.out.println("Couldn't introspect " + c.getName());
                System.exit(1);
            }

        	PropertyDescriptor[] properties = bi.getPropertyDescriptors();
            int nextRootIndex = 1;

    		// initialize property boxes
    		for(int i = 0; i < numItems; i++) {

    			String myPrefix = boxes[i].getPrefix();
    			String myTitle  = boxes[i].getTitle();

    			PropertyDescriptor myProp = null;

                // Pull property specific to this box from array
        		for(int prop = 0; prop < properties.length;  prop++) {
    	    		String propertyName = properties[prop].getName();

    	    		if (propertyName.startsWith(myPrefix) &&
    	    			propertyName.substring(myPrefix.length()).indexOf(myTitle) != -1) {

    	    		    myProp = properties[prop];
    	        		boxes[i].setProperty(myProp);
    	        		break;
    	    		}
    	    	}

                // Compile root index
                int j;
                for(j = 0; j < i;  j++) {
                    PropertyDescriptor property = boxes[j].getProperty();
                    if (myPrefix.equals(boxes[j].getPrefix())) {
                        String subString = property.getName().substring(myPrefix.length(), property.getName().length());
                        if (subString.indexOf(myTitle) != -1) {
                            rootIndex[i] = rootIndex[j];
                            break;
                        }
                    }
                }

                if(j == i && i != 0) {
                    rootIndex[i] = nextRootIndex++;
                }

    		}
        }

    	private int getIndex(String title) {
    		int index = -1;

    		for(int i = 0; i < numItems;  i++) {
    			if (title.equals(boxes[i].getTitle())) {
    				index = i;
    				break;
    			}
    		}
    		return index;
    	}

    	private ModifierGeneral getModifier(String title) {
    		ModifierGeneral mod = null;

    		int index = getIndex(title);
    		if (index != -1) {
    			mod = boxes[index].getModifier();
    		}
    		return mod;
    	}

    	public boolean hasRead(String title) {
    		boolean value = false;
    		
            int index = getIndex(title);
    		if (index != -1) {
    			value = boxes[index].hasRead();
    		}
    		// Probably should throw an exception if title does
    		// not match the title of any of the boxes
    		return value;
    	}

    	public boolean hasWrite(String title) {
    		boolean value = false;
    		
            int index = getIndex(title);
    		if (index != -1) {
    		    value = boxes[index].hasWrite();
    		}
    		// Probably should throw an exception if title does
    		// not match the title of any of the boxes
    		return value;
    	}

    	public int getRootIndex(String title) {
    		int index = -1;
            int idx = getIndex(title);
    			if (idx != -1) {
    				index = rootIndex[idx];
    		}
            return index;
    	}

    } // End LatticeEditorBoxPropertyArray

    /**
     * LatticeEditorBoxProperty
     * @author rrassler
     *
     */
    private class LatticeEditorBoxProperty {

    	private String prefix    = null;
    	private String title     = null;
    	private boolean hasRead  = false;
    	private boolean hasWrite = false;
        private PropertyDescriptor property = null;
    	private ModifierGeneral modifier = null;

    	public LatticeEditorBoxProperty(String title, String prefix, PropertyDescriptor prop) {
    		this.title   = title;
    		this.prefix  = prefix;

    		propertyUpdate(prop);
    	}

    	private void propertyUpdate(PropertyDescriptor prop) {

    		this.property = prop;

    		if (property == null) {
    			hasRead = false;
    			hasWrite = false;
    		}
    		else {

				modifier = new ModifierGeneral(currentLattice.getPrimitive(),
						                       property.getName());
                modifier.setLabel(this.title);

        		if (property.getReadMethod() == null)
    			    hasRead = false;
    	        else
    	        	hasRead = true;

    		    if (property.getWriteMethod() == null)
			        hasWrite = false;
	            else
	    	        hasWrite = true;   			
    		}

    	}

    	public void setProperty(PropertyDescriptor prop) {
    		propertyUpdate(prop);
    	}

    	public boolean hasRead() {
    		return hasRead;
    	}

    	public boolean hasWrite() {
    		return hasWrite;
    	}

    	public String getTitle() {
    		return title;
    	}

    	public ModifierGeneral getModifier() {
    		return modifier;
    	}

    	public String getPrefix() {
    		return prefix;
    	}
    	
        public PropertyDescriptor getProperty() {
	        return property;
        }

    } // End LatticeEditorBoxProperty

}
