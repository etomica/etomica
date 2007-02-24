package etomica.modules.crystalviewer;
import java.awt.Component;
import java.awt.FlowLayout;
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

import etomica.action.Action;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DeviceBox;
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

/**
 * Class that produces a panel with controls that permit editing of
 * the structure of a lattice.  Permits changing of type primitive and
 * parameters (sizes and angles) for lattice vectors.
 */
 
 /* History
  * 09/07/02 (DAK) new
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
                updateBoxes();
                update();
            }
        });

        anglePanel = new JPanel(new GridLayout(1,3));
        sizePanel = new JPanel(new GridLayout(1,3));
        
        DeviceSlider nSlider = new DeviceSlider(null, new NModifier());
        nSlider.setPrecision(0);
        nSlider.setMinimum(1);
        nSlider.setMaximum(10);
        nSlider.getSlider().setValue(DEFAULT_SIZE);
        nSlider.setLabel("Size");
        nSlider.setShowBorder(true);
        
        boxPanel = new JPanel(new GridLayout(2,1));
        boxPanel.add(sizePanel);
        boxPanel.setBorder(new javax.swing.border.TitledBorder("Primitive vectors"));
       
        updateBoxes();
        
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
    
    protected void updateBoxes() {
        for (int i=0; i<angleBoxes.length; i++) {
            anglePanel.remove(angleBoxes[i].graphic());
        }
        angleBoxes = new DeviceBox[0];
        for (int i=0; i<sizeBoxes.length; i++) {
            sizePanel.remove(sizeBoxes[i].graphic());
        }
        if (angleBoxes.length > 0) {
            boxPanel.remove(anglePanel);
        }
        sizeBoxes = new DeviceBox[0];
        
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
        boolean hasAngle = false;
        for (int i=0; i<properties.length; i++) {
            String name = properties[i].getName();
            if (!name.startsWith("size") && !name.startsWith("angle") && !name.equals("cubicSize") || properties[i].getWriteMethod() == null) {
                // we don't actually want setSize as it wants an array
                continue;
            }
            ModifierGeneral modifier = new ModifierGeneral(primitive, name);
            DeviceBox newBox = new DeviceBox();
            if (name.startsWith("angle")) {
                // angle
                modifier.setLabel(name.substring(5));
                hasAngle = true;
                newBox.setUnit(Degree.UNIT);
                anglePanel.add(newBox.graphic());
                angleBoxes = (DeviceBox[])Arrays.addObject(angleBoxes, newBox);
            }
            else {
                // might be A, B, C, AB, CubicSize
                if (name.startsWith("size")) {
                    modifier.setLabel(name.substring(4));
                }
                sizePanel.add(newBox.graphic());
                sizeBoxes = (DeviceBox[])Arrays.addObject(sizeBoxes, newBox);
            }
            newBox.setModifier(modifier);
            newBox.setPostAction(new Action() {
                public void actionPerformed() {
                    update();
                }
                public String getLabel() {
                    // because they deserve it.
                    return null;
                }
            });
        }
        if (hasAngle) {
            boxPanel.setLayout(new GridLayout(2,1));
            boxPanel.add(anglePanel);
        }
        else {
            boxPanel.setLayout(new FlowLayout());
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
    }

}