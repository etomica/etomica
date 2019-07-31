/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.clustergenerator;

import java.awt.BorderLayout;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;

import etomica.virial.cluster.ClusterDiagram;
import etomica.virial.cluster.ClusterGenerator;

/**
 * @author andrew
 * @author modified by skkwak
 *
 */

public class MyApplet extends javax.swing.JApplet {
    public JPanel myPanel;
    protected ClusterDiagram cluster;
    protected ClusterGenerator generator;
    protected int nBody;
    protected TextArea output;
    protected TextField nBodyField, nRootField;
    protected Checkbox allPermutationsCB;
    protected Checkbox onlyDoublyConnectedCB, onlySinglyConnectedCB;
    protected Checkbox excludeNodalPointCB;
    protected Checkbox excludeArticulationPointCB, excludeArticulationPairCB;
    protected Checkbox generateGraphicsCB;
    protected Checkbox eDrawBondsCB, fDrawBondsCB;
    protected Checkbox reeHooverCB;
    protected boolean running = false;
    protected JButton goButton;
    public Starter starter = new Starter();
    protected Stopper stopper = new Stopper();
    protected GraphicsToggler graphicsToggler = new GraphicsToggler();
    protected boolean stopNow;
    protected JScrollPane graphicOut;
    protected JPanel graphicOut2;

    private JPanel graphicPanel;
    protected JPanel managerPanel;
    protected JTabbedPane displayPanel;
    protected GridBagConstraints gbConst;    
    protected TextArea textAreaForBonds;
    protected GridBagLayout gbLayout;
    protected Checkbox checkBox, checkBox1, checkBox2;
    protected JRadioButton radioB, radioB1;
    protected BondTypeToggler bondTypeToggler;
	protected JTextField  totalNCluster;
    private boolean firstTimeDrawn = true;
    protected boolean allPermutations;
    
    public void init() {
        Dimension d = this.getSize();
        d.width = 650; d.height = 670;
        setSize(d);
    	myPanel = new JPanel(); // the basis panel of applet
        myPanel.setSize(d);
        myPanel.setLayout(new BorderLayout());
        myPanel.setBorder(new TitledBorder(
                new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
                    ,"  Cluster Diagram Generator  ",TitledBorder.CENTER,TitledBorder.TOP
                    ,new Font(null,Font.BOLD,16),java.awt.Color.blue));	
        myPanel.addComponentListener(new ComponentEventMyPanel());
        getContentPane().add(myPanel);

    	gbLayout = new GridBagLayout();
    	gbConst = new GridBagConstraints();

    	managerPanel = new JPanel();// the basis of the manager panel
	    managerPanel.setLayout(gbLayout);
	    managerPanel.setBorder(new TitledBorder(
	            new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
	                ," Cluster Diagram Manager ",TitledBorder.LEFT,TitledBorder.TOP
	                ,new Font(null,Font.BOLD,14),Color.black));

	    goButton = new JButton("Start"); // start button is created
	    goButton.setBackground(Color.lightGray);
	    d.width = 80; d.height = 50;
	    goButton.setPreferredSize(d);
	    goButton.addActionListener(starter);
//	    gbConst.weightx = 0; gbConst.weighty = 0;
//        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(goButton, managerPanel, 0,0,1,1); // start button is added to the manager panel
        
        JPanel pointSelectorPanel = new JPanel(new GridLayout(2,2)); // panel of number of points is created 
        pointSelectorPanel.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
                ,"Number of Points",TitledBorder.LEFT,TitledBorder.TOP
                ,new Font(null,Font.BOLD,12),Color.black));
        pointSelectorPanel.add(new Label("TOTAL",Label.CENTER)); 
        pointSelectorPanel.add(new Label("ROOT", Label.CENTER));
        nBodyField = new TextField("", 2);
        pointSelectorPanel.add(nBodyField);
        nRootField = new TextField("0", 2);
        pointSelectorPanel.add(nRootField);
	    gbConst.weightx = 2; gbConst.weighty = 0;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(pointSelectorPanel, managerPanel, 0,1,1,1); // the panel of number of points is added to the manager panel       

        JPanel clusterTypePanel = new JPanel(new GridLayout(3,1)); // option panel of graphics of strings
        clusterTypePanel.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
                ,"Display Option ",TitledBorder.LEFT,TitledBorder.TOP
                ,new Font(null,Font.BOLD,12),Color.black));
        checkBox = new Checkbox("Show Points");
        checkBox1 = new Checkbox("Show Number : White Points");
        checkBox2 = new Checkbox("Show Number : Black Points");
        checkBox.setState(true);
        checkBox1.setState(false);checkBox2.setState(false);
        checkBox.addItemListener(new PointsToggler());
        checkBox1.addItemListener(new NumbersToggler());
        checkBox2.addItemListener(new NumbersToggler());
        clusterTypePanel.add(checkBox);
        clusterTypePanel.add(checkBox1);
        clusterTypePanel.add(checkBox2);
        gbConst.weightx = 1; gbConst.weighty = 0;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(clusterTypePanel, managerPanel, 0,2,1,1); // clusterTypePanel is added to the manager panel        

        JTabbedPane bondType = new JTabbedPane(); // any bond types can be added in this pane
        JPanel bondTypeSelectorPanel = new JPanel(new GridLayout(1,2)); // options for array of bonds
        bondTypeSelectorPanel.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
                ,"",TitledBorder.LEFT,TitledBorder.TOP
                ,new Font(null,Font.BOLD,12),Color.black));
        bondTypeToggler = new BondTypeToggler();
        radioB = new JRadioButton("Default");
        radioB.setSelected(true);
        radioB1 = new JRadioButton("Java");
        radioB.setSelected(false);
        radioB.addItemListener(bondTypeToggler);
        radioB1.addItemListener(bondTypeToggler);
        bondTypeSelectorPanel.add(radioB);
        bondTypeSelectorPanel.add(radioB1);
        JPanel bondInteractionType = new JPanel(new GridLayout(1,3)); // options for interaction of bonds (etc. e-bonds)
        bondInteractionType.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
                ,"",TitledBorder.LEFT,TitledBorder.TOP
                ,new Font(null,Font.BOLD,12),Color.black));
        
        fDrawBondsCB = new Checkbox("Draw f-bonds");
        fDrawBondsCB.setState(true);
        eDrawBondsCB = new Checkbox("Draw e-bonds");
        fDrawBondsCB.addItemListener(new BondInteractionTypeToggler());
        eDrawBondsCB.addItemListener(new BondInteractionTypeToggler());
        bondInteractionType.setToolTipText("Must enable Ree-Hoover diagrams");
        bondInteractionType.add(fDrawBondsCB);
        bondInteractionType.add(eDrawBondsCB);
        fDrawBondsCB.setEnabled(false);
        eDrawBondsCB.setEnabled(false);
        bondType.add("Bond Type", bondInteractionType);
        bondType.add("Bond Array", bondTypeSelectorPanel); // bondTypeSelectorPanel is added to bondType panel
	    gbConst.weightx = 2; gbConst.weighty = 0;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(bondType, managerPanel, 0,3,1,1);
        
        JPanel temp = new JPanel(new GridLayout(4,2)); // any constraints for cluster generation can be added here
        temp.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.red, Color.blue) 
                ,"Constraints",TitledBorder.LEFT,TitledBorder.TOP
                ,new Font(null,Font.BOLD,12),Color.black));
        allPermutationsCB = new Checkbox("All Permutations");
        allPermutationsCB.setBackground(Color.lightGray);
        onlySinglyConnectedCB = new Checkbox("Singly Connected Clusters");        
        onlySinglyConnectedCB.setBackground(Color.lightGray);
        onlyDoublyConnectedCB = new Checkbox("Doubly Connected Clusters");        
        onlyDoublyConnectedCB.setBackground(Color.lightGray);
        excludeNodalPointCB = new Checkbox("Exclude strings with a nodal point");
        excludeNodalPointCB.setBackground(Color.lightGray);
        excludeArticulationPointCB = new Checkbox("Exclude strings with an articulation point");
        excludeArticulationPointCB.setBackground(Color.lightGray);
        excludeArticulationPairCB = new Checkbox("Exclude strings with an articulation pair");
        excludeArticulationPairCB.setBackground(Color.lightGray);
        reeHooverCB = new Checkbox("Use Ree-Hoover diagrams");
        reeHooverCB.setBackground(Color.lightGray);
        reeHooverCB.setEnabled(false);
        ReeHooverToggler reeHooverToggler = new ReeHooverToggler();
        reeHooverCB.addItemListener(reeHooverToggler);
        ConnectedToggler connectedToggler = new ConnectedToggler();
        onlySinglyConnectedCB.addItemListener(connectedToggler);
        onlyDoublyConnectedCB.addItemListener(connectedToggler);
        onlyDoublyConnectedCB.addItemListener(reeHooverToggler);
/*        generateGraphicsCB = new Checkbox("Generate Graphics");
        generateGraphicsCB.setBackground(Color.lightGray);
        generateGraphicsCB.setState(true);
        generateGraphicsCB.addItemListener(graphicsToggler);
        temp.add(generateGraphicsCB);*/
        temp.add(allPermutationsCB);
        temp.add(excludeNodalPointCB);
        temp.add(onlySinglyConnectedCB);
        temp.add(excludeArticulationPointCB);
        temp.add(onlyDoublyConnectedCB);
        temp.add(excludeArticulationPairCB);
        temp.add(reeHooverCB);
	    gbConst.weightx = 0; gbConst.weighty = 0;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(temp, managerPanel, 1,0,4,1); // constraints panel is added to the manager panel       

        displayPanel = new JTabbedPane(); // displayPanel is created to add cluster graphics and bonds array to myPanel
//        displayPanel.setFocusable(true);
        makeGraphicsPanel(); // every components related to graphics of cluster and bonds are created and arranged 
        
        output = new TextArea("", 20, 80); // textArea for the array of bonds
        output.setBackground(Color.white);
        output.setEditable(false);
        displayPanel.add("Cluster Diagram Bonds", output); // added to displayPanel

        myPanel.add(managerPanel, BorderLayout.SOUTH);// managerPanel is added to myPanel
        myPanel.add(displayPanel, BorderLayout.CENTER);// displayPanel is added to myPanel
        
    }

    
    public class Starter implements ActionListener {
        private int totalNumberOfCluster;
        private boolean drawPoints = true;
        private boolean drawNumbersWhite = false;
        private boolean drawNumbersBlack = false;
        
        public synchronized void actionPerformed(ActionEvent evt) {
            running = true;
            stopNow = false;
            goButton.setText("Stop");
            goButton.removeActionListener(starter);
            goButton.addActionListener(stopper);
            totalNumberOfCluster = 0;
            totalNCluster.setText("");
            textAreaForBonds.setText("Bonds of Selected Cluster :");
            textAreaForBonds.append("\n\n");
            textAreaForBonds.append("Weight of Selected Cluster : ");            
            Thread runner = new Thread(new Runnable() {
                public void run() {foo();}
            });
            runner.start();
        }
        
        public void foo() {
            String nBodyText = nBodyField.getText();
            int nRootPoints;
            String nRootText = nRootField.getText();
            try {
                nBody = Integer.valueOf(nBodyText).intValue();
            }
            catch (NumberFormatException e) {
                nBodyField.requestFocus();
                allDone();
                return;
            }
            try {
                nRootPoints = Integer.valueOf(nRootText).intValue();
            }
            catch (NumberFormatException e) {
                nRootField.requestFocus();
                allDone();
                return;
            }
            
            cluster = new ClusterDiagram(nBody,nRootPoints);
            generator = new ClusterGenerator(cluster);
            allPermutations = allPermutationsCB.getState();
            generator.setAllPermutations(allPermutations);
            generator.setOnlyConnected(onlySinglyConnectedCB.getState());
            generator.setOnlyDoublyConnected(onlyDoublyConnectedCB.getState());
            generator.setExcludeArticulationPoint(excludeArticulationPointCB.getState());
            generator.setExcludeArticulationPair(excludeArticulationPairCB.getState());
            generator.setExcludeNodalPoint(excludeNodalPointCB.getState());
            generator.setMakeReeHover(onlyDoublyConnectedCB.getState() && reeHooverCB.getState());
            cluster.reset();
            generator.reset();
            if (reeHooverCB.isEnabled() && reeHooverCB.getState()) {
                generator.calcReeHoover();
            }
            output.setText("");
            graphicOut2.removeAll();
            output.append(""+cluster.mReeHooverFactor);
            if (!allPermutations) {
                output.append("/"+cluster.mNumIdenticalPermutations);
            }
            output.append("x\t");
            output.append(cluster.toString() + "\n");
            if (graphicOut != null) {
                addCluster(cluster);
            } else { totalNumberOfCluster+=1;
                totalNCluster.setText(String.valueOf(totalNumberOfCluster));
            }
           
            while (!stopNow && generator.advance()) {
                output.append(""+cluster.mReeHooverFactor);
                if (!allPermutations) {
                    output.append("/"+cluster.mNumIdenticalPermutations);
                }
                output.append("x\t");
                output.append(cluster.toString() + "\n");
                if (graphicOut != null) {
//                    graphicOut2.setLayout(new GridLayout(0,2));
                    addCluster(cluster);
                } else { totalNumberOfCluster+=1;
                         totalNCluster.setText(String.valueOf(totalNumberOfCluster));
                }
            }
            allDone();
        }
        
        public synchronized void allDone() {
            running = false;
            stopNow = false;
            goButton.setText("Start");
            goButton.removeActionListener(stopper);
            goButton.addActionListener(starter);
            myPanel.updateUI();
//            System.out.println(graphicOut2.getWidth()+" "+graphicOut2.getHeight());            
        }

        public void addCluster(ClusterDiagram newCluster) {
        	totalNumberOfCluster+=1;
        	totalNCluster.setText(String.valueOf(totalNumberOfCluster));
            ClusterDiagram copy = new ClusterDiagram(newCluster);
            boolean fBonds = true;
            boolean eBonds = false;
            if (fDrawBondsCB.isEnabled()) {
                fBonds = fDrawBondsCB.getState();
                eBonds = eDrawBondsCB.getState();
            }
            ClusterPanel cp = new ClusterPanel(copy, drawPoints, drawNumbersWhite, drawNumbersBlack, fBonds, eBonds);
            //cp.addMouseListener(new MouseEventClusterPanel());
            Dimension d = new Dimension();
            d.width=100;
            d.height=100;
            cp.setSize(d);
            JPanel jp = new JPanel(new java.awt.GridLayout(1,1));
            cp.addMouseListener(new MouseEventClusterPanel());
            jp.setPreferredSize(d);
            jp.setMinimumSize(d);
            jp.add(cp);
            jp.setBackground(null);
            graphicOut2.add(jp);
            myPanel.updateUI();
        }
        public void setDrawPoints(boolean b){
        	drawPoints = b;
        }
        public void setDrawNumbersOfWhite(boolean b){
        	drawNumbersWhite = b;
        }
        public void setDrawNumbersOfBlack(boolean b){
        	drawNumbersBlack = b;
        }

    }
    
    public class Stopper implements ActionListener {
        public synchronized void actionPerformed(ActionEvent e) {
            if (running) {
                stopNow = true;
                goButton.setText("Waiting...");
                goButton.removeActionListener(stopper);
            }
        }
    }
    
    public class ConnectedToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
            if (evt.getStateChange() == ItemEvent.SELECTED) {
                if (evt.getSource() == onlyDoublyConnectedCB) {
                    onlySinglyConnectedCB.setState(false);
                }
                else {
                    onlyDoublyConnectedCB.setState(false);
                }
            }
        }
    }
    
    public class GraphicsToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
            if (evt.getStateChange() == ItemEvent.SELECTED) {
                makeGraphicsPanel();
            }
            else {
//            	myPanel.remove(graphicOut);
                graphicOut = null;
                myPanel.updateUI();
            }
        }
    }

    public void makeGraphicsPanel() {
    	graphicPanel = new JPanel();
    	graphicPanel.setLayout(gbLayout);
    	
    	graphicOut2 = new JPanel();
    	graphicOut2.setLayout(new GridLayout(0,(int)((myPanel.getWidth()-60)/100.0),0,0));
        graphicOut2.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.black, Color.gray)));	   
        graphicOut = new JScrollPane();
//        graphicOut2.addMouseListener(new MouseEventClusterPanel());
        graphicOut.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.black, Color.gray)));	   

        graphicOut2.setPreferredSize(null);
        graphicOut.setBackground(null);
        graphicOut2.setBackground(null);
        graphicOut.getViewport().setPreferredSize(null);
        graphicOut.getViewport().setBackground(Color.white);
        graphicOut.getViewport().add(graphicOut2);
        gbConst.weightx = 1; gbConst.weighty = 1;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(graphicOut,graphicPanel, 0,0,4,1);

        JPanel tmp = new JPanel();
        textAreaForBonds = new TextArea("Bonds of Selected Cluster :", 4, 65);
        textAreaForBonds.append("\n\n");
        textAreaForBonds.append("Weight of Selected Cluster : ");
        textAreaForBonds.setBackground(Color.white);
        textAreaForBonds.setEditable(false);
        tmp.add(textAreaForBonds);
        gbConst.weightx = 0; gbConst.weighty = 0;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(tmp,graphicPanel, 1,0,3,1);

        JPanel textPanel = new JPanel(new GridLayout(2,0));
        textPanel.setBorder(new TitledBorder(new EtchedBorder(EtchedBorder.RAISED, Color.black, Color.gray)));
        JTextField  jtf = new JTextField();
        jtf.setBackground(Color.lightGray);
        jtf.setText("Total Number of Cluster");
        jtf.setEditable(false);
        textPanel.add(jtf);
        totalNCluster = new JTextField();
        textPanel.add(totalNCluster);
        gbConst.weightx = 0; gbConst.weighty = 0;
        gbConst.fill = GridBagConstraints.BOTH;
        addComponent(textPanel,graphicPanel, 1,3,1,1);
        
        if(firstTimeDrawn) {
        	displayPanel.add("Cluster Diagram Graphics", graphicPanel); firstTimeDrawn= false;        
        } else {
        	displayPanel.remove(0);displayPanel.add(graphicPanel, "Cluster Diagram Graphics",  0);
        }
    }

    /**
     * Componet c is added to JPanel panel following the rule of GridBagLayout manager
     */
    public void addComponent(Component c, JPanel panel, int row, int column, int width, int height) {
    	gbConst.gridx = column;gbConst.gridy = row;
    	gbConst.gridwidth = width;gbConst.gridheight = height;
    	gbLayout.setConstraints(c, gbConst);
    	panel.add(c);
    }
    /**
     * This class is for showing bonds of selected cluster in GraphicPanel
     * @author skkwak
     */
	public class MouseEventClusterPanel implements MouseListener {
		public void mouseClicked(MouseEvent arg0) {
            ClusterDiagram thisCluster = ((ClusterPanel)arg0.getComponent()).getCluster();
            textAreaForBonds.setText("Bonds of Selected Cluster : \n"
            		+thisCluster.toString()+"\n");
            textAreaForBonds.append("Weight of Selected Cluster : "+thisCluster.mReeHooverFactor);
            if(!allPermutations){
                textAreaForBonds.append(" / "+thisCluster.mNumIdenticalPermutations);
            }
		}
		
		public void mouseEntered(MouseEvent arg0) {}
		public void mouseExited(MouseEvent arg0) {}
		public void mousePressed(MouseEvent arg0) {}
		public void mouseReleased(MouseEvent arg0) {}
	}
	/**
	 * This class is to be notified to graphicOut2 when the window of applet is resized
	 * @author skkwak
	 */
    public class ComponentEventMyPanel implements ComponentListener{
    	public void componentHidden(ComponentEvent e){}
    	public void componentMoved(ComponentEvent e){}
    	public void componentShown(ComponentEvent e){}
    	public void componentResized(ComponentEvent e){
            int n = (int)((graphicOut2.getWidth()-50)/100.0); // 100 is the dimension of a drawn cluster
            graphicOut2.setLayout(new GridLayout(0,n));
            graphicOut2.updateUI();
    	}    	        
    }
    /**
     * Option of showing particles or not
     * @author skkwak
     */
    public class PointsToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
            if (evt.getStateChange() == ItemEvent.SELECTED) {
                starter.setDrawPoints(true);
            }
            else {
            	starter.setDrawPoints(false);
            }
        }
    }
    /**
     * Options for showing index of particle or not
     * @author skkwak
     */
    public class NumbersToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
        	if(evt.getSource().equals(checkBox1)){
	            if (evt.getStateChange() == ItemEvent.SELECTED) {
	                starter.setDrawNumbersOfWhite(true);
	            }
	            else {
	            	starter.setDrawNumbersOfWhite(false);
	            }
        	} else if(evt.getSource().equals(checkBox2)){
	            if (evt.getStateChange() == ItemEvent.SELECTED) {
	                starter.setDrawNumbersOfBlack(true);
	            }
	            else {
	            	starter.setDrawNumbersOfBlack(false);
	            }
        	}
        }
    }
    /**
     * Options for type of bond array
     * @author skkwak
     *
     * TODO Type of bond can be created in GenCluster class or here. It's up to you, Andrew :)
     */
    public class BondTypeToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
        	if(evt.getSource().equals(radioB)){
				if (evt.getStateChange() == ItemEvent.SELECTED) {
					radioB1.setSelected(false);
					radioB1.updateUI();
		            textAreaForBonds.setText(" Not Yet!! :)");
				}
        	}
        	if(evt.getSource().equals(radioB1)){
				if (evt.getStateChange() == ItemEvent.SELECTED) {
					radioB.setSelected(false);
					radioB.updateUI();
		            textAreaForBonds.setText(" Not Yet!! :)");
				}
        	}
        	
        }
    }
    /**
     * Options for type of bond-interaction  (e-bonds, f-bonds, etc)
     * @author skkwak
     *
     * TODO we do it later!!!
     */
    public class BondInteractionTypeToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
        	if(evt.getSource().equals(fDrawBondsCB)){
				if (evt.getStateChange() != ItemEvent.SELECTED) {
                    eDrawBondsCB.setState(true);
				} 
        	}
            else if (evt.getSource().equals(eDrawBondsCB)){
                if (evt.getStateChange() != ItemEvent.SELECTED) {
                    fDrawBondsCB.setState(true);
                }
            }
        }
    }

    /**
     * Options for type of bond-interaction  (e-bonds, f-bonds, etc)
     * @author skkwak
     *
     * TODO we do it later!!!
     */
    public class ReeHooverToggler implements ItemListener {
        public void itemStateChanged(ItemEvent evt) {
            boolean enabled = onlyDoublyConnectedCB.getState();
            reeHooverCB.setEnabled(enabled);
            enabled = enabled && reeHooverCB.getState();
            eDrawBondsCB.setEnabled(enabled);
            fDrawBondsCB.setEnabled(enabled);
        }
    }
}

