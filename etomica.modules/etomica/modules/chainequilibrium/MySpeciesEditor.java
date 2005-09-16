/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import etomica.action.Action;
import etomica.action.ActionGroup;
import etomica.atom.Atom;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
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
		
		// ********** Marker ********
		System.out.println("Made a New Species Editor");
		
		//nSlider.setDisplayPhase(DisplayPhase1);
        int majorSpacing = 50;
		nSlider.setMinimum(0);
		nSlider.setMaximum(100);
        nSlider.getSlider().setLabelTable(nSlider.getSlider().createStandardLabels(majorSpacing));
        
        
//        nSlider.addChangeListener(new javax.swing.event.ChangeListener() {
//			
//			public void stateChanged(javax.swing.event.ChangeEvent evt) {
//
//				// ********* Marker 
//				//System.out.println("Graphic: State Change");
//				
//				AtomIteratorListSimple iter = new AtomIteratorListSimple(
//						sim.phase1.speciesMaster.atomList);
//				iter.reset();
//				while (iter.hasNext()) {
//					//            			System.out.println(iter.peek().toString());
//                    Atom[] a = (Atom[])iter.nextAtom().allatomAgents[sim.idx];
//                    a[0] = null;
//                    a[1] = null;
//				}
//				sim.integratorHard1.reset();
//                simGraphic.displayPhase1.repaint();
//			}
//		});
		//            nSlider = new DeviceSlider(new NMoleculeModulator(s));
		//            nSlider.setShowBorder(true);
		//// nSlider.setLabel(label);
		//			nSlider.setLabel("Atom count");
		//	        nSlider.setMinimum(0);
		//	        nSlider.setMaximum(40);
		//	        nSlider.getSlider().setSnapToTicks(true);
		//	        nSlider.graphic(null).setSize(new java.awt.Dimension(40,30));

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
				((AtomTypeLeaf)species.type.getSpecies().getFactory().getType()).setMass(newMass);
				try {
                    sim.integratorHard1.reset();
                } catch (ConfigurationOverlapException e) {
                }
//				 ********* Marker 
				System.out.println("Graphic: Action Preformed");
			}
		};
		
//		 ********* Marker 
		System.out.println("Graphic: Added Mass Listener");
		
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
                    AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.phase1);
                    iter.reset();
                    while (iter.hasNext()) {
                        //                      System.out.println(iter.peek().toString());
                        Atom[] a = (Atom[])iter.nextAtom().allatomAgents[sim.idx];
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
            targetAction = new ActionGroup(new Action[] {targetAction, anotherAction});
        }
    }
    
} //end of MySpeciesEditor

