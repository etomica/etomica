/**
 * DefineConfigurationFrame
 *
 * The DefineConfigurationFrame class is responsible for creating a new JInternalFrame that gives the user
 * the ability to determine the proper Configuration of atoms in the the AtomTypeArray to make the newly
 * created, user defined Molecule.
 *
 * @author Bryan C. Mihalick
 * 9/18/00
 */
 
package simulate.gui;

import simulate.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.DefaultComboBoxModel;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JInternalFrame;
import javax.swing.JPanel;
import javax.swing.JTextField;


public class DefineConfigurationFrame extends JInternalFrame {
    JInternalFrame mainFrameHandle;
    String configuration;
    String molTitle;
    static int molNameCounter = 0;
    JButton periodicTable = new JButton("Periodic Table");
    JButton createMolecule = new JButton("Create Molecule");
    JComboBox configList = new JComboBox();
    JPanel mainPanel = new JPanel();
    JPanel configPanel = new JPanel();
    JTextField molName = new JTextField(8);
    AtomType[] typeArray;
    SimulationEditor simulationEditor;
    
    DefineConfigurationFrame(SimulationEditor ed){
        simulationEditor = ed;
        mainFrameHandle = this;
        setResizable(true);
        setMaximizable(true);
        setIconifiable(true);
        setClosable(false);
        setVisible(true);
        setBounds(515, 60, 275, 100);
        setTitle("Configure New Molecule");
        
        DefaultComboBoxModel atomTypeList = new DefaultComboBoxModel();
        configList.setModel(atomTypeList);
        configList.addItem(new String("Sequential"));
        configList.addItem(new String("Tetrahedral"));
        configList.addItem(new String("Ring"));
        configList.setSelectedIndex(0);
        configList.addItemListener(new ItemListener(){
            public void itemStateChanged(ItemEvent evt) {
	            configuration = ((String)((JComboBox)evt.getSource()).getSelectedItem());
                
                if (configuration == "Sequential")
	                System.out.println("Sequential");
	            else if (configuration == "Tetrahedral")
	                System.out.println("Tetrahedral");
	            else if (configuration == "Ring")
	                System.out.println("Ring");
            }});
        mainPanel.add(configList);
        
        molName.setText("Untitled" + molNameCounter);
        mainPanel.add(molName);
        
        
        createMolecule.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e){
                molTitle = molName.getText();
                if (molTitle.startsWith("Untitled"))
                    molNameCounter++;
                try {
                    mainFrameHandle.setClosed(true);
                }
                catch (java.beans.PropertyVetoException pve) {}
                if (typeArray[0] instanceof AtomType.Disk) {
                    SpeciesUserDefinedDisks newSpeciesDisk = new SpeciesUserDefinedDisks(typeArray);
                    newSpeciesDisk.setName(molTitle);
                    simulationEditor.speciesEditor.componentList.addElement(newSpeciesDisk);
                }
                else { 
                    AtomType.Wall[] wallArray = new AtomType.Wall[1];
                    for (int i = 0; i < typeArray.length; i++)
                        wallArray[i] = (AtomType.Wall)typeArray[i];
                    SpeciesWalls newSpeciesWall = new SpeciesWalls(1,(wallArray));
                    newSpeciesWall.setName(molTitle);
                    simulationEditor.speciesEditor.componentList.addElement(newSpeciesWall);
                }
                simulationEditor.potential1Editor.setEnabled(true);
                simulationEditor.potential2Editor.setEnabled(true);
                simulationEditor.integratorEditor.setEnabled(true);
                simulationEditor.phaseEditor.setEnabled(true);
                simulationEditor.controllerEditor.setEnabled(true);
                simulationEditor.displayEditor.setEnabled(true);
                simulationEditor.deviceEditor.setEnabled(true);
                simulationEditor.meterEditor.setEnabled(true);
                
                //add to Species Tabbed Pane RadioButton list    
            }});
        mainPanel.add(createMolecule);
        
        periodicTable.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e){
                Etomica.DesktopFrame.desktop.add(new PeriodicTable());    
            }});
        mainPanel.add(periodicTable);
        
        getContentPane().add(mainPanel);
    }// end of DefineConfigurationFrame constructor
    
    public void setAtomTypes(AtomType[] a) { typeArray = ((AtomType[])a.clone()); }
    public AtomType[] getAtomTypes() { return typeArray; }
}//end of DefineConfigurationFrame class