package etomica.graphics;
import etomica.*;
import etomica.units.*;

//Java2 imports
//import java.util.Iterator;

import etomica.utility.java2.Iterator;

/**
 * Permits selection of temperature from a discrete set of values.  Also has option
 * to select adiabatic.  Enforces selection by manipulating properties of integrator.
 * By default, last integrator added to simulation is the only one controlled by this device.
 *
 * @author David Kofke
 */
public class DeviceThermoSelector extends Device implements EtomicaElement {
    
    public String getVersion() {return "DeviceThermoSelector:01.03.23/"+Device.VERSION;}
    /**
     * Descriptive text label to be displayed with the value
     */
    private javax.swing.JComboBox selector;
    private javax.swing.JLabel label;
    private final String adiabaticString = "Adiabatic";
    private Integrator integrator;
    private javax.swing.JPanel panel;
    
    private boolean includeAdiabatic = true;
    
    public DeviceThermoSelector() {
        this(Simulation.instance);
    }
    public DeviceThermoSelector(Simulation sim) {
        super(sim);
        selector = new javax.swing.JComboBox(new Object[] {new Object()}); //swingall JComboBox doesn't chokes if we don't give initialize without an object in the list
        setTemperatures(new double[] {200.0, 400.0, 600.0});
        selector.setEditable(false);
        label = new javax.swing.JLabel("");
        setUnit(new PrefixedUnit(Kelvin.UNIT));
        
        panel = new javax.swing.JPanel(new java.awt.BorderLayout(0,1));
        panel.add(label, java.awt.BorderLayout.NORTH);
        panel.add(selector, java.awt.BorderLayout.SOUTH);
        panel.setBorder(new javax.swing.border.EmptyBorder(3,3,3,3));
        
        
/*        selector.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent e){
                Object item = ((javax.swing.JComboBox)e.getSource()).getSelectedItem();
		        if(item==adiabaticString) integrator.setIsothermal(false);
		        else {
		            integrator.setIsothermal(true);
		            integrator.setTemperature(unit.toSim(((Double)item).doubleValue()));
		        }
            }
        });*/
        selector.addItemListener(new java.awt.event.ItemListener() {
		    public void itemStateChanged(java.awt.event.ItemEvent event) {
		        if(integrator==null) return;
		        Object item = event.getItem();
		        if(item==adiabaticString) integrator.setIsothermal(false);
		        else {
		        	integrator.pause();
		            integrator.setIsothermal(true);
		            integrator.setTemperature(unit.toSim(((Double)item).doubleValue()));
		            integrator.reset();
		            integrator.unPause();
		        }
		    }
		});//end of addItemListener
		
		//add mediator so that last integrator added to simulation is controlled by this device
        sim.elementCoordinator.addMediatorPair(new MediatorGraphic.DeviceIntegrator(sim.elementCoordinator) {
            public void add(Integrator integrator) {
                DeviceThermoSelector.this.integrator = integrator;
            }
            public void add(Device device) {
                if(device != DeviceThermoSelector.this) return;
                for(Iterator ip=mediator.parentSimulation().getIntegratorList().iterator(); ip.hasNext(); ) {
                    Integrator integrator = (Integrator)ip.next();
                    if(integrator.wasAdded())  {//will make last integrator the one
                        DeviceThermoSelector.this.integrator = integrator;
                    }
                }
            }
        });
		
    }//end of constructor
    
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
        selector.setEnabled(true);
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
}