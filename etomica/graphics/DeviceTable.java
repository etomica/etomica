package etomica.graphics;
import etomica.*;
import etomica.units.*;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.JPanel;
import javax.swing.Box;
import javax.swing.JScrollPane;

/**
 * Presents a table of numeric properties that can be edited by typing in values.
 *
 * @author David Kofke
 */
 
 //does not work work Etomica because isn't called with constructor that 
 //provides modulators.
 //need to return a non-null graphic (instantiate panel in all constructors)
 //and set up set/get methods for modulators
public class DeviceTable extends Device /*implements EtomicaElement*/ {
    
    public String getVersion() {return "DeviceTable:01.04.17/"+Device.VERSION;}

    public JTable table;
    MyTableData dataSource;
    Modulator[] modulators;
    DimensionedDoubleEditor[] editors;
    PropertyText[] views;
    JPanel panel;

    public DeviceTable() {
        this(Simulation.instance);
    }
    
    public DeviceTable(Simulation sim) {
        super(sim);
    }
    
    public DeviceTable(Simulation sim, Modulator[] mods) {
        super(sim);
        modulators = mods;
        editors = new DimensionedDoubleEditor[modulators.length];
        views = new PropertyText[modulators.length];
        setupTable();
        for(int i=0; i<modulators.length; i++) {
            editors[i] = new DimensionedDoubleEditor(modulators[i].getDimension());
            editors[i].setValue(modulators[i].getValue());
            new Adapter(editors[i], modulators[i]);
            views[i] = new PropertyText(editors[i]);
        }
    }
    
    private void setupTable() {
        panel = new JPanel();
        panel.setSize(100,150);
        dataSource = new MyTableData();   //inner class, defined below
        table = new JTable(dataSource);
        panel.add(new JScrollPane(table));
        table.getColumn("Units").setCellEditor(new UnitEditor());
        table.getColumn("Value").setCellEditor(new ValueEditor());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Editable table of property values");
        return info;
    }

    public void setPropertyLabels(String[] labels) {
        int n = labels.length;
        if(n != modulators.length) {
            System.out.println("Warning:  Number of labels given to DeviceTable disagrees with number of modulators");
            n = Math.min(n, modulators.length);
        }
        for(int i=0; i<n; i++) {
            modulators[i].setLabel(labels[i]);
        }
    }
    
    public java.awt.Component graphic(Object obj) {return panel;}

    class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        
        MyTableData() {
            columnNames = new String[] {"Property", "Value", "Units"};
            columnClasses = new Class[] {String.class, String.class, Unit.class};
        }
        
        public Object getValueAt(int row, int column) {
            switch(column) {
                case 0: return modulators[row].getLabel();
                case 1: return editors[row].getAsText();
                case 2: return editors[row].getUnit();
                default: return null;
            }
        }
        
        public int getRowCount() {return modulators.length;}
        public int getColumnCount() {return 3;}
       
        public boolean isCellEditable(int row, int column) {return column != 0;}        
        public String getColumnName(int column) {return columnNames[column];}
        public Class getColumnClass(int column) {return columnClasses[column];}
    }

    /**
     * ComboBox editor of units column of table.
     */
    class UnitEditor extends javax.swing.DefaultCellEditor {
        
        public UnitEditor() {
            super(new javax.swing.JTextField());  //dummy argument since there is no default constructor in superclass
            setClickCountToStart(1);
        }
        public java.awt.Component getTableCellEditorComponent(
                JTable table, Object value, boolean isSelected, int row, int column) {
            return editors[row].unitSelector();
        }
        public boolean isCellEditable(java.util.EventObject e) {return true;}
        
    }//end of UnitEditor
    
    /**
     * TextBox editor of values column of table.
     */
    class ValueEditor extends javax.swing.DefaultCellEditor {
        
        public ValueEditor() {
            super(new javax.swing.JTextField());  //dummy argument since there is no default constructor in superclass
            setClickCountToStart(1);
        }
        public java.awt.Component getTableCellEditorComponent(
                JTable table, Object value, boolean isSelected, int row, int column) {
            return views[row];
        }
        public boolean isCellEditable(java.util.EventObject e) {return true;}
        
    }//end of ValueEditor
    
    /**
     * Class to tie together DimensionedDoubleEditor, modulator, and table.
     */
    class Adapter implements java.beans.PropertyChangeListener {
        private DimensionedDoubleEditor editor;
        private Modulator modulator;
        Adapter(DimensionedDoubleEditor ed, Modulator mod) {
            editor = ed;
            modulator = mod;
            ed.addPropertyChangeListener(this);
        }
        //method called when DimensionedDoubleEditor changes value
        public void propertyChange(java.beans.PropertyChangeEvent evt) {
            double newValue = ((Double)editor.getValue()).doubleValue();
            modulator.setValue(newValue);
            panel.repaint();
        }        
    }//end of Adapter
    
    /**
     * main method to demonstrate and test this class.
     */
/*    public static void main(String[] args) {
        
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);

        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;

        PotentialSquareWell potential = new PotentialSquareWell(sim);
        sim.p2 = new P2SimpleWrapper(sim,sim.potential);
        Modulator mod1 = new Modulator(sim.integrator, "timeStep");
        Modulator mod2 = new Modulator(potential, "epsilon");
        DeviceTable table = new DeviceTable(Simulation.instance, new Modulator[] {mod1, mod2});
        sim.integrator.setIsothermal(true);
        sim.integrator.setTemperature(Kelvin.UNIT.toSim(300.));
		DisplayBox box1 = new DisplayBox();
		DisplayBox box2 = new DisplayBox();
		box1.setDatumSource(mod1);
		box2.setDatumSource(mod2);
		
		Simulation.instance.elementCoordinator.go(); 
		                                    
        f.getContentPane().add(Simulation.instance.panel());         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
 */   
}

