package etomica.modules.crystalviewer;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.util.HashMap;

import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.Action;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceBoxValueChangedEvent;
import etomica.graphics.DeviceBoxValueChangedListener;
import etomica.graphics.DeviceSlider;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Primitive;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.species.Species;
import etomica.units.Degree;
import etomica.units.Quantity;
import etomica.util.Arrays;
import etomica.graphics.DeviceBoxValueChangedListener;

/**
 * Class that produces a panel with controls that permit editing of
 * the structure of a lattice.  Permits changing of type primitive and
 * parameters (sizes and angles) for lattice vectors.
 */
public class LatticeEditor {
    
    public static int DEFAULT_SIZE = 6;
    protected int size;
    protected BravaisLattice currentLattice;
    private DeviceBox[] angleBoxes, sizeBoxes;
    public JPanel anglePanel, sizePanel;
    public JPanel boxPanel;
    protected CrystalViewer viewer;
    protected final HashMap latticeNameHash;
    protected PrimitiveVectorBox pvBox;
	private final String[] fieldTitles = {"A", "B", "C",
                                          "Alpha", "Beta", "Gamma" };
    private final String[] fieldPrefix = {"size", "size", "size",
                                          "angle", "angle", "angle"};
    private static final int UNIT_CELL_START_INDEX = 0;
    private static final int UNIT_CELL_END_INDEX = 2;
    private static final int ANGLE_START_INDEX = 3;
    private static final int ANGLE_END_INDEX = 5;

    public LatticeEditor(CrystalViewer viewer, BravaisLattice[] lattices, String[] latticeNames) {
    
        this.viewer = viewer;
        phase = viewer.phase;
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
                initBoxes();
                update();
            }
        });

        anglePanel = new JPanel(new GridLayout(1,3));
        sizePanel = new JPanel(new GridLayout(1,3));
        sizePanel.setBorder(new TitledBorder("Cell Dimensions"));
        DeviceSlider nSlider = new DeviceSlider(null, new NModifier());
        nSlider.setPrecision(0);
        nSlider.setMinimum(1);
        nSlider.setMaximum(10);
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

        // Create a primitive vector box instance
        pvBox = new PrimitiveVectorBox(this.fieldTitles, this.fieldPrefix);

    	// Unit cell size entry fields
        for(int box = UNIT_CELL_START_INDEX; box <= UNIT_CELL_END_INDEX; box++) {

            DeviceBox newBox = new DeviceBox();

            sizePanel.add(newBox.graphic());
            sizeBoxes = (DeviceBox[])Arrays.addObject(sizeBoxes, newBox);

            newBox.setPostAction(new Action() {
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

            newBox.setPostAction(new Action() {
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

        int numAtoms = size*size*size;
        if (currentLattice instanceof BravaisLatticeCrystal) {
            numAtoms *= ((BravaisLatticeCrystal)currentLattice).getBasis().getScaledCoordinates().length;
        }
        IVector dimensions = phase.getBoundary().getDimensions();
        dimensions.E(currentLattice.getPrimitive().getSize());
        dimensions.TE(size);
        phase.setDimensions(dimensions);
        phase.getAgent(species).setNMolecules(numAtoms);
        ConfigurationLattice config = new ConfigurationLattice(currentLattice);
        config.initializeCoordinates(phase);
        viewer.update(currentLattice);
    }

    public JPanel getPanel() {return panel;}

    private JPanel panel;
    protected Phase phase;
    protected Species species;

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
        public etomica.units.Dimension getDimension() {return Quantity.DIMENSION;}
        public String getLabel() {
            return "a label";
        }
    } // end protected class NModifier



    private class PrimitiveVectorBox {

        private class ModifierCreation {
        	private String[] titles = null;
        	private String[] prefix = null;
        	private boolean[] hasRead = null;
        	private boolean[] hasWrite = null;
        	private ModifierGeneral[] modifiers = null;
            private int numItems = 0;
          
        	ModifierCreation(String[] title, String[] prefix,
        			PropertyDescriptor[] properties, Primitive primitive) {

        	    if (title.length == prefix.length) {
        	    	this.numItems = title.length;
            		this.titles = title;
            		this.prefix = prefix;
            		this.hasRead = new boolean[numItems];
            		this.hasWrite = new boolean[numItems];
            		for (int i = 0; i < numItems; i++) {
                		this.hasRead[i] = false;
                		this.hasWrite[i] = false;
            		}

            		modifiers = new ModifierGeneral[this.numItems];
        	    }
        	    else {
        	    	return;
        	    }

        	    for(int i = 0; i < this.numItems;  i++) {

    	    		for(int prop = 0; prop < properties.length;  prop++) {
        	    		String propertyName = properties[prop].getName();

        	    		if (propertyName.startsWith(this.prefix[i])) {
        	    			if (propertyName.substring(this.prefix[i].length()).indexOf(this.titles[i]) != -1) {

        	    				modifiers[i] = new ModifierGeneral(primitive, propertyName);

        	                    modifiers[i].setLabel(this.titles[i]);
                                if (properties[prop].getReadMethod() == null)
                                	hasRead[i] = false;
                                else
                                	hasRead[i] = true;
                                if (properties[prop].getWriteMethod() == null)
                                	hasWrite[i] = false;
                                else
                                	hasWrite[i] = true;
                                break;
        	    			}
        	    		}
        	    	}
        	    }
        	}

        	private int getIndex(String title) {
        		int index = -1;
        	    for (int i = 0; i < this.numItems;  i++) {
        	    	if (title.equals(this.titles[i])) {
        	    		index = i;
        	    		break;
        	    	}
        	    }
        	    return index;
        	}

        	public ModifierGeneral getModifier(String fieldTitle) {
        		
        	    ModifierGeneral mod = null;
        	    int index = getIndex(fieldTitle);

        	    if (index != -1) {
        	        mod = modifiers[index];
        	    }
        		return mod;
        	}
        	
        	public boolean getHasRead(String fieldTitle) {
        		boolean read = false;
        	    int index = getIndex(fieldTitle);

        	    if (index != -1) {
        	    	read = hasRead[index];
        	    }
        		return read;
        	}

        	public boolean getHasWrite(String fieldTitle) {
        		boolean write = false;
        	    int index = getIndex(fieldTitle);

        	    if (index != -1) {
        	    	write = hasWrite[index];
        	    }
        		return write;
        	}
        
        } // private class ModifierCreation



    	private PropertyDescriptor[] properties = null;  
    	private ModifierCreation modCreater = null;
    	private String[] fieldTitles;
        private String[] fieldPrefix;

    	PrimitiveVectorBox(String[] titles, String[] prefix) {
    		this.fieldTitles = titles;
    		this.fieldPrefix = prefix;
            initialize();
    	}

    	protected void initialize() {
    		properties = null;
    		modCreater = null;
            
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
            properties = bi.getPropertyDescriptors();

            modCreater = new ModifierCreation(fieldTitles, fieldPrefix,
        			properties, primitive);
    	}

    	private ModifierGeneral getModifier(String title) {
    		return modCreater.getModifier(title);
    	}

    	private boolean hasRead(String title) {
    		return modCreater.getHasRead(title);
    	}

    	private boolean hasWrite(String title) {
    		return modCreater.getHasWrite(title);
    	}
    } // private class PrimitiveVectorBox

}