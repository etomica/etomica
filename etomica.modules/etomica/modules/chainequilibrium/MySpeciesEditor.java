/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceNSelector;
import etomica.species.Species;

/**
 * @author Matt Moynihan
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

//=================================================================
//panel containing species-editing devices
class MySpeciesEditor extends javax.swing.JPanel {

	//	public DeviceSlider nSlider;
	public DeviceNSelector nSlider;
	public Species species;
	boolean initializing;
    ChainEquilibriumSim sim;
	public final javax.swing.JTextField mass = new javax.swing.JTextField("40");

	
	public MySpeciesEditor(final ChainEquilibriumGraphic simGraphic, Box box, Species s, String label) {
		super();
		species = s;
        sim = simGraphic.sim;
		nSlider = new MyNSelector(simGraphic, box, species);
		
		//nSlider.setDisplayBox(DisplayBox1);
        int majorSpacing = 50;
		nSlider.setMinimum(0);
		nSlider.setMaximum(100);
        nSlider.getSlider().setLabelTable(nSlider.getSlider().createStandardLabels(majorSpacing));
        
		//listener for changes to mass textbox
		java.awt.event.ActionListener myListener = new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent event) {
				if (initializing)
					return;
				int value;
				try {
					value = Integer.parseInt(mass.getText());
				} catch (NumberFormatException ex) {
					return;
				}
				if (value < 1)
					value = 1;
				if (value > 1000000)
					value = 1000000;
				final double newMass = value;
				mass.setText(Integer.toString(value));
				((ElementSimple)((AtomTypeLeaf)species.getMoleculeType()).getElement()).setMass(newMass);
				try {
                    sim.integratorHard1.reset();
                } catch (ConfigurationOverlapException e) {
                }
			}
		};
		
		mass.addActionListener(myListener);
		mass.setBorder(new javax.swing.border.TitledBorder("Mass"));
		mass.setColumns(6);
		mass.setOpaque(false);
		setLayout(new java.awt.FlowLayout());
		add(nSlider.graphic(null));
		add(mass);
		setBorder(new javax.swing.border.TitledBorder(label));
	}
    
    class MyNSelector extends DeviceNSelector {
        MyNSelector(final ChainEquilibriumGraphic simGraphic, Box box, Species species) {
            super(simGraphic.sim.getController());
            setResetAction(new SimulationRestart(simGraphic.sim));
            setBox(box);
            setSpecies(species);
            
            Action anotherAction = new Action() {
                public void actionPerformed() {
                    AtomAgentManager agentManager = sim.getAgentManager();
                    AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.box);
                    iter.reset();
                    for (IAtom atom = iter.nextAtom(); atom != null;
                         atom = iter.nextAtom()) {
                        //                      System.out.println(iter.peek().toString());
                        IAtom[] a = (IAtom[])agentManager.getAgent(atom);
                        for(int i = 0; i < a.length; i++) {
                            a[i] = null;
                        }
                    }
                    try {
                        sim.integratorHard1.reset();
                    } catch (ConfigurationOverlapException e) {
                    }
                    simGraphic.getDisplayBox(sim.box).repaint();
                    
                }
                
            };
            targetAction = new ActionGroupSeries(new Action[] {targetAction, anotherAction});
        }
    }
    
} //end of MySpeciesEditor

