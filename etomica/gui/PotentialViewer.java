package simulate.gui;

import simulate.Potential;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import javax.swing.DefaultListModel;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.JButton;
import javax.swing.JInternalFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

public class PotentialViewer extends JInternalFrame {
    static DefaultListModel listModel = new DefaultListModel();
    static JList potentialList = new JList(listModel);
    static JPanel panel = new JPanel();
    static JScrollPane scroller = new JScrollPane(potentialList);    
    public static Potential[][] defPotArray;
    public static String title;
    public static JButton ok = new JButton("OK");
    
    /**
     * Holds all constraints needed for displaying the next awt or swing component
     */
    final GridBagConstraints gbc = new GridBagConstraints();
    
    /**
     * Determines how to display an awt or swing component by using gbc from above
     */
    final GridBagLayout gbl = new GridBagLayout();
    
    /**
     * Envelopes a simulation component in order to have the component's properties listed in the 
     * property sheet.
     */
    public static Wrapper wrapper = null;
    
    /**
     * Internal frame that lists the properties of a component
     */
    public static PropertySheet propSheet;
    
    public SimulationEditor simulationEditor;
    
    PotentialViewer(SimulationEditor ed){
        simulationEditor = ed;
        setBounds(515, 60, 250, 200);
        setResizable(true);
        setVisible(false);
        panel.setLayout(gbl);
        
        potentialList.addListSelectionListener(new javax.swing.event.ListSelectionListener(){
            public void valueChanged(javax.swing.event.ListSelectionEvent lse){
                wrapper = new Wrapper(potentialList.getSelectedValue(), title, "simulate.gui." + title); 
                propSheet.setTarget(wrapper);
	            propSheet.addInternalFrameListener(new InternalFrameAdapter(){
	                public void internalFrameClosed( InternalFrameEvent ife ){
	                    potentialList.clearSelection();
	                }});
                propSheet.setVisible(true);
                propSheet.repaint();
            }});
		potentialList.getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        potentialList.setSize(250,140);
        gbc.fill = gbc.BOTH;
        gbc.weightx = 95;
        gbc.weighty = 95;
        gbl.addLayoutComponent(scroller, gbc);
        
        ok.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent evt){
                Etomica.DesktopFrame.desktop.getComponent(0).setVisible(false);
	            simulationEditor.potential1Editor.rightPaneList.clearSelection();
	            simulationEditor.potential2Editor.rightPaneList.clearSelection();
            }});
        gbc.weightx = 0;
        gbc.weighty = 0;
        gbc.gridx++;
        gbc.fill = gbc.NONE;
        gbc.anchor = gbc.CENTER;
        gbl.addLayoutComponent(ok, gbc);
        
    	scroller.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        panel.add(scroller);
        panel.add(ok);
        getContentPane().add(panel);
        Etomica.DesktopFrame.desktop.add(this);
    }//End of PotentialViewer Constructor
    
    public static void setParam(Potential[][] potArray, PropertySheet p, String t) {
        defPotArray = potArray;
        propSheet = p;
        title = t;
        listModel.removeAllElements();
        for (int i = 0; i < defPotArray.length; i++) {
            for (int j = 0; j < defPotArray[i].length; j++) {
                listModel.addElement(defPotArray[i][j]);
            }
        }
        
    }
}//End of PotentialViewer Class