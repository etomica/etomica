package etomica;
import java.io.*;
import etomica.units.*;
import javax.swing.*;
import java.lang.reflect.*;
import java.beans.*;
//Java2 imports
//import java.util.Iterator;

import etomica.utility.Iterator;


/**
 * Implements logging functionality, performing writing of simulation information and results to a log file.
 *
 * @author David Kofke
 */

public class DisplayToLog extends Display {
    
    public String getVersion() {return "DisplayToLog:01.03.02.0/"+Display.VERSION;}
    
    private PrintWriter logFile;
    
    public DisplayToLog() {
        this(Simulation.instance);
    }
    public DisplayToLog(Simulation sim) {
        this(sim, "default.log");
    }
    public DisplayToLog(String fileName) {
        this(Simulation.instance, fileName);
    }
    public DisplayToLog(Simulation sim, String fileName) {
        super(sim);
        try {
            logFile = new PrintWriter(
                        new BufferedWriter(
                            new FileWriter(fileName)));
        }
        catch(IOException ioe) {
            System.err.println("I/O Exception in DisplayToLog");
    	    System.err.println("Cannot open file: " + fileName);
        }
    }
    
    public java.awt.Component graphic(Object obj) {return null;}
    
    public void intervalAction(Integrator.IntervalEvent evt) {
        if(evt.type() == Integrator.IntervalEvent.START) writeHeading();
        else if(evt.type() == Integrator.IntervalEvent.INTERVAL && --iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
	    }
	    else if(evt.type() == Integrator.IntervalEvent.DONE) writeEnding();
    }
    
    private void writeHeading() {
        logFile.println("Heading");
        for(Iterator iter=parentSimulation().allElements().iterator(); iter.hasNext(); ) {
            printProperties(iter.next());
        }
        logFile.println("=========================");
    }
    
    private void writeEnding() {
        logFile.println("Ending");
        logFile.close();//should do this when controller announces end
    }
    
    public void doUpdate() {
        for(Iterator iter=parentSimulation().meterList().iterator(); iter.hasNext(); ) {
            MeterAbstract meter = (MeterAbstract)iter.next();
            if(!(meter instanceof Meter)) continue;
            logFile.print(((Meter)meter).average() + " " + ((Meter)meter).error()+ "  ");
        }
        logFile.println();
    }
    
    public void printProperties(Object object) {
        
        PropertyDescriptor[] properties = null;

        //Introspection to get array of all properties
        BeanInfo bi = null;
        try {
	        bi = Introspector.getBeanInfo(object.getClass());
	        properties = bi.getPropertyDescriptors();
	    } 
	    catch (IntrospectionException ex) {
	        System.out.println("DisplayLog: Couldn't introspect " + ex);
	        return;
	    }

	    logFile.println("*******************");
	    logFile.println(object.toString());

        //Loop through all properties and print values
	    for (int i = 0; i < properties.length; i++) {
	        // Don't display hidden or expert properties.
	        if (properties[i].isHidden() || properties[i].isExpert())
		        continue;

	        Object value = null;
	        JLabel label = null;
	        Dimension dimension = null;

	        String name = properties[i].getDisplayName();  //Localized display name 
	        Class type = properties[i].getPropertyType();  //Type (class) of this property
	        Method getter = properties[i].getReadMethod(); //method used to read value of property in this object
	        // Only display read/write properties.
	        if (getter == null) continue;
	        if (type == Class.class) continue;
	        if (Dimension.class.isAssignableFrom(type)) continue;
	        if (name.equals("name")) continue;

	        try {
	            //read the current value of the property
		        Object args[] = { };
		        try {value = getter.invoke(object, args);}
		        catch(NullPointerException ex) {value = null;}

        	        //property is a dimensioned number
        	        if(type == Double.TYPE) {
        	            //try to get dimension from get(property)Dimension() method
        	            dimension = etomica.units.Dimension.introspect(object,name,bi);
        	            //try to get dimension from getDimension() method
                        if(dimension == null) dimension = etomica.units.Dimension.introspect(object,"",bi);
        	        }
	            
	        } //end of try
	        catch (InvocationTargetException ex) {
		        System.err.println("Skipping property " + name + " ; exception on target: " + ex.getTargetException());
		        ex.getTargetException().printStackTrace();
		        continue;
	        } 
	        catch (Exception ex) {
		        System.err.println("Skipping property " + name + " ; exception: " + ex);
		        ex.printStackTrace();
		        continue;
	        }
	        if(value == null) continue;
	        Unit unit;
	        if(type == Double.TYPE && dimension != null) {
    	        unit = dimension.defaultIOUnit();
    	        double v = ((Double)value).doubleValue();
    	        logFile.println(name + " " + unit.fromSim(v) + " " + unit.toString());
	        }
	        else logFile.println(name + " " + value.toString());
	    }//end of loop over properties
    }//end of printProperties
    
    //tests and demonstrates use of this class
 /*   public static void main(String[] args) {
        
        Default.ATOM_SIZE = 1.0;
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        sim.setUnitSystem(new UnitSystem.LJ());
        
        Species species = new SpeciesDisks(32);
        Potential potential = new PotentialHardDisk();
        Potential2 p2 = new P2SimpleWrapper(potential);
        Integrator integrator = new IntegratorHard();
        Controller controller = new Controller();
        Meter meter = new MeterPressureHard();
        DisplayBox box = new DisplayBox();
        Phase phase = new Phase();
        DisplayPhase display = new DisplayPhase();
        
        DisplayToLog displayLog = new DisplayToLog();
                
        box.setMeter(meter);
        box.setPrecision(8);
        
        sim.elementCoordinator.go();
        phase.setDensity(0.5);
        display.setScale(2.0);
        controller.setMaxSteps(100);

        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        f.add(sim);         //access the static instance of the simulation to
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
    */
}