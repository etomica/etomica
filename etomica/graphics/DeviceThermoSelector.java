package etomica.graphics;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.Action;
import etomica.Controller;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.simulations.HSMD3D;
import etomica.units.PrefixedUnit;

/**
 * Permits selection of temperature from a discrete set of values.  Also has option
 * to select adiabatic.  Enforces selection by manipulating properties of integrator.
 * By default, last integrator added to simulation is the only one controlled by this device.
 *
 * @author David Kofke
 */
public class DeviceThermoSelector extends Device implements EtomicaElement {
    
     public DeviceThermoSelector(Controller controller, final Integrator integrator) {
        this();
        setController(controller);
        setIntegrator(integrator);
     }
     
     public DeviceThermoSelector() {
        super();
        selector = new javax.swing.JComboBox(new Object[] {new Object()});
        setTemperatures(new double[] {200.0, 400.0, 600.0});
        selector.setEditable(false);
        label = new javax.swing.JLabel("");
        setUnit(Default.UNIT_SYSTEM.temperature());

        panel = new javax.swing.JPanel(new java.awt.BorderLayout(0,1));
        panel.add(label, java.awt.BorderLayout.NORTH);
        panel.add(selector, java.awt.BorderLayout.SOUTH);
        panel.setBorder(new javax.swing.border.EmptyBorder(3,3,3,3));

        
        //listener to combo box gets value and initiates action
        selector.addItemListener( new ItemListener() {
            public void itemStateChanged(ItemEvent event) {
                Object item = event.getItem();
                if(item==adiabaticString) {
                    isothermal = false;
                } else {
                    isothermal = true;
                    temperature = unit.toSim(((Double)item).doubleValue());
                }
                doAction(targetAction);
            }
        });
       
    }//end of constructor
    
    public void setIntegrator(final Integrator integrator) {
        this.integrator = integrator;
        targetAction = new Action() {
            public void actionPerformed() {
                if(integrator == null) return;
                integrator.setIsothermal(isothermal);
                if(isothermal) {
                    integrator.setTemperature(temperature);
                }
                integrator.reset();
            }
            public String getLabel() {
                return DeviceThermoSelector.this.getLabel().toString();
            }
        };
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Select isothermal at a set of temperatures, or adiabatic");
        return info;
    }
    
    public javax.swing.JComboBox getSelector() {return selector;}
    public javax.swing.JLabel getLabel() {return label;}
    
    public final void setIncludeAdiabatic(boolean b) {
        if(selector.getItemCount() == 0) includeAdiabatic = false;
        if(includeAdiabatic && !b) selector.removeItem(adiabaticString);//remove if was there and removing
        else if(!includeAdiabatic && b) selector.addItem(adiabaticString);//add if wasn't there and adding
        includeAdiabatic = b;
    }
    public boolean isIncludeAdiabatic() {return includeAdiabatic;}
    
    /**
     * Sets the i-th item (counting from 0) in the list as the one currently selected.
     * Argument outside bounds is interpreted as first or last item on list.
     */
    public void setSelected(int i) {
        if(i < 0) i = 0;
        if(i >= selector.getItemCount()) i = selector.getItemCount()-1;
        selector.setSelectedIndex(i);
    }
    
    public void setTemperatures(double[] t) {
        selector.removeAllItems();
        setIncludeAdiabatic(includeAdiabatic);
        for(int i=0; i<t.length; i++) {
            selector.addItem(new Double(t[i]));
        }
    }

    /**
     * Override superclass setUnit method to update label when unit is changed
     */
    public void setUnit(PrefixedUnit u) {
        super.setUnit(u);
        label.setText("Temperature ("+unit.symbol()+")");
    }
    
    
    /**
     * Returns the GUI element for display in the simulation.
     * Consists of a combo box used for the selector.
     */
    public java.awt.Component graphic(Object obj) {
        return panel;
       // return selector;
    }
    
    private ItemListener itemListener;
    private javax.swing.JComboBox selector;
    private javax.swing.JLabel label;
    private final String adiabaticString = "Adiabatic";
    private javax.swing.JPanel panel;
    
    private boolean includeAdiabatic = true;
    private Integrator integrator;
    private boolean isothermal = false;
    private double temperature;
    private Action targetAction;
    
    //main method to demonstrate and test class
    public static void main(String[] args) {
        
        HSMD3D sim = new HSMD3D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        
        DeviceThermoSelector device = new DeviceThermoSelector(sim.getController(), sim.integrator);
        device.setTemperatures(new double[] {0.5, 1.0, 2.0, 5.0});
        graphic.add(device);
        LinkedList displays = graphic.displayList();
        for(Iterator iter=displays.iterator(); iter.hasNext();) {
            Display next = (Display)iter.next();
            if(next instanceof DisplayPhase) {
                ((DisplayPhase)next).setColorScheme(new ColorSchemeTemperature(0.5, 5.0));
            }
        }
        graphic.makeAndDisplayFrame();
    }//end of main

}