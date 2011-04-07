package etomica.virial.GUI.containers;

import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;

import javax.swing.border.SoftBevelBorder;



public class ParameterTabs extends JPanel{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public JTabbedPane Tabs;
	public Species TabPanel1;
	private RunParam TabPanel2;
	private JPanel TabPanel3;
	private JPanel TabPanel4;
	private String LabelText;
	
	
	public ParameterTabs(int NoOfTabs, int Width){
		super(new GridLayout(1, 1));
		
		Tabs = new JTabbedPane();
		Tabs.setTabPlacement(SwingConstants.LEFT);
		Tabs.setAlignmentX((float)0.5);
		Tabs.setAlignmentY((float)0.5);
		Tabs.setPreferredSize(new Dimension(437,616));
		Tabs.setBorder(new SoftBevelBorder(SoftBevelBorder.RAISED));
		//MakeIndividualTabs(this.Tabs,NoOfTabs);
	}
	
	public void MakeIndividualTabs(JTabbedPane Tabs, int NoOfTabs, Species species, RunParam runparam){
		
		for(int Count = 0; Count < NoOfTabs;Count++){
			
			if(Count == 0){
				LabelText = "Species A";
				//TabPanel1 = new Species();
				this.setTabPanel1(species);
				Tabs.addTab(LabelText, this.getTabPanel1());
				
			}
			if(Count == 1){
				LabelText = "Species B";
				//TabPanel2 = new RunParam();
				this.setTabPanel2(runparam);
				Tabs.addTab(LabelText, this.getTabPanel2());
			}
			if(Count == 2){
				LabelText = "Run Parameters";
			this.TabPanel3 = new JPanel();
			Tabs.addTab(LabelText,TabPanel3);
				
			}
			
			if(Count == 2){
				LabelText = "Moves";
			this.TabPanel4 = new JPanel();
			Tabs.addTab(LabelText,TabPanel4);
				
			}
			add(Tabs);
		}
	}
	
	public JTabbedPane getTabs() {
		return Tabs;
	}

	public void setTabs(JTabbedPane tabs) {
		Tabs = tabs;
	}

	public Species getTabPanel1() {
		return TabPanel1;
	}

	public void setTabPanel1(Species tabPanel1) {
		TabPanel1 = tabPanel1;
	}

	public RunParam getTabPanel2() {
		return TabPanel2;
	}

	public void setTabPanel2(RunParam tabPanel2) {
		TabPanel2 = tabPanel2;
	}

	public String getLabelText() {
		return LabelText;
	}

	public void setLabelText(String labelText) {
		LabelText = labelText;
	}
	
	

}


