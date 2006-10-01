/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.atom.Atom;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.chem.elements.ElementSimple;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceNSelector;

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
	public SpeciesAgent species;
	boolean initializing;
    ReactionEquilibrium sim;
	public final javax.swing.JTextField mass = new javax.swing.JTextField("40");

	
	public MySpeciesEditor(final ReactionEquilibriumGraphic simGraphic, SpeciesAgent s, String label) {
		super();
		species = s;
        sim = simGraphic.simulation;
		nSlider = new MyNSelector(simGraphic, species);
		
		//nSlider.setDisplayPhase(DisplayPhase1);
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
				((ElementSimple)((AtomTypeLeaf)species.type.getSpecies().getMoleculeType()).getElement()).setMass(newMass);
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
        MyNSelector(final ReactionEquilibriumGraphic simGraphic, SpeciesAgent species) {
            super(simGraphic.simulation, species);
            
            Action anotherAction = new Action() {
                public void actionPerformed() {
                    Atom[][] agents = sim.getAgents(sim.phase1);
                    AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.phase1);
                    iter.reset();
                    while (iter.hasNext()) {
                        //                      System.out.println(iter.peek().toString());
                        Atom[] a = agents[iter.nextAtom().getGlobalIndex()];
                        a[0] = null;
                        a[1] = null;
                    }
                    try {
                        sim.integratorHard1.reset();
                    } catch (ConfigurationOverlapException e) {
                    }
                    simGraphic.displayPhase1.repaint();
                    
                }
                public String getLabel() {return "";}
                
            };
            targetAction = new ActionGroupSeries(new Action[] {targetAction, anotherAction});
        }
    }
    
} //end of MySpeciesEditor

