package etomica;
import etomica.units.*;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.JPanel;
import javax.swing.Box;
import javax.swing.JScrollPane;

/**
 * Presents a table of numeric properties that can be edited by typing in values.
 */
 
public class DeviceTable extends Device implements EtomicaElement {
    
    public JTable table;
    MyTableData dataSource;
    Modulator[] modulators;
    DimensionedDoubleEditor[] editors;
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
        setupTable();
        for(int i=0; i<modulators.length; i++) {
            editors[i] = new DimensionedDoubleEditor(modulators[i].getDimension());
            editors[i].setValue(modulators[i].getValue());
            new Adapter(editors[i], modulators[i]);
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

    public void repaint() {table.repaint();}
    
    class MyTableData extends AbstractTableModel {
        
        String[] columnNames;
        Class[] columnClasses;
        
        MyTableData() {
            columnNames = new String[] {"Property", "Value", "Units"};
            columnClasses = new Class[] {String.class, Double.class, Unit.class};
        }
        
        public Object getValueAt(int row, int column) {
            switch(column) {
                case 0: return modulators[row].getLabel();
                case 1: return new Double(modulators[row].getValue());
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
            return editors[row].valueEditor();
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
 //           double newValue = ((Double)((DimensionedDoubleEditor)evt.getSource()).getValue()).doubleValue();
            double newValue = ((Double)editor.getValue()).doubleValue();
            modulator.setValue(newValue);
            table.repaint();
        }        
    }//end of Adapter
    
    /**
     * main method to demonstrate and test this class.
     */
    public static void main(String[] args) {
        
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);

	    IntegratorHard integratorHard1 = new IntegratorHard();
	    SpeciesDisks speciesDisks1 = new SpeciesDisks();
	    Phase phase1 = new Phase();
	    Potential potential = new PotentialSquareWell();
	    P2SimpleWrapper p2 = new P2SimpleWrapper(potential);
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();

        Modulator mod1 = new Modulator(integratorHard1, "timeStep");
        Modulator mod2 = new Modulator(potential, "epsilon");
        DeviceTable table = new DeviceTable(Simulation.instance, new Modulator[] {mod1, mod2});
                                            
		Simulation.instance.elementCoordinator.go(); 
		                                    
        f.add(Simulation.instance.panel());         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}

