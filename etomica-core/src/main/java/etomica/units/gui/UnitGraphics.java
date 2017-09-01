/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.gui;

import etomica.units.CompoundUnit;
import etomica.units.Unit;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Force;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

//import java.lang.Number;

public class UnitGraphics {

	// The target dimension.
	private Dimension targetDimension;

	// The final unit.
	private Unit currentCompoundUnit;

	// The button that can be pressed when the current compound unit has the
	// same dimension as the target dimension.
	private JButton submitChoice;

	// Label that displays the current action.
	private JLabel label;

	// Label that displays the target dimension.
	private JLabel targetDim;

	// A button for adding units to the used units list.
	private JButton addUnitButton;

	// Button for removing units from the used units list.
	private JButton removeUnitButton;

	// Button to increase the exponent that the selected unit has.
	private JButton raisePower;

	// Button that lowers the exponent for the selected unit.
	private JButton lowerPower;

	// A combo box for selecting what dimension units the units list contains.
	private JComboBox dimensionChoice;

	// A combo box that allows selection of a unit system to filter the unit
	// list.
	private JComboBox unitSystemChoice;

	// The list of the units that are currently being used.
	private JList usedUnits;

	// The list off all the units.
	private JList allUnits;

	// The center panel that contains most of the buttons and the combo boxes.
	private JPanel centerPane;

	// The top panel that contains the label.
	private JPanel topPane;

	// The bottom panel that contains the submit unit button.
	private JPanel bottomPane;

	// The right list that has a set of units based on the filters used.
	private JScrollPane rightListPane;

	// The left list that displays the selected units and the power each is
	// raised to.
	private JScrollPane leftListPane;

	// A string that contains whatever the selected unit is on the used units
	// list.
	private String selectedUsedUnit;

	// A string that contains the selected unit on the possible units list.
	private String selectedUnit;

	// A string array containing all the dimensions.
	private String[] allDimensionStrings = Lister.dimensionString();

	// A string array containing all the units.
	private String[] allUnitStrings = Lister.unitString();

	// A string array containing all unit systems.
	private String[] allUnitSystemStrings = Lister.unitSystemString();

	// The amount of rows for the list panels to display.
	private int rowCount = 15;

	// A vector containing all the units that are currently used.
	private Vector unitsVector = new Vector();

	// A vector with all the exponents for the used units.
	private Vector exponentVector = new Vector();

	// The list that is displayed for the used units.
	private DefaultListModel usedUnitsModel = new DefaultListModel();

	// The list that is displayed for the chosen set of units.
	private DefaultListModel selectedUnitsSet = new DefaultListModel();

	private String dl;

	private int di;

	private boolean switching = false;

	// The object used to check actions from buttons and combo boxes.
	private ActionCheck acheck = new ActionCheck();

	/*
	 * The main function creates and shows the GUI.
	 */
	public void main(String[] args) {
		targetDimension = Force.DIMENSION;
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				createAndShowGUI();
			}
		});

	}

	public Unit startWithDim(Dimension dim) {
		targetDimension = dim;

		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				createAndShowGUI();
			}
		});

		return currentCompoundUnit;
	}

	public void setDim(Dimension dim) {
		targetDimension = dim;
	}

	/*
	 * A function that creates a label and the panel it goes into.
	 */
	private Component createTop() {
		label = new JLabel("Choose Your Unit");
		targetDim = new JLabel("Target Dimension:");
		targetDim.setText("Target Dimension: " + targetDimension.toString());
		topPane = new JPanel(new GridLayout(0, 1));
		topPane.add(label);
		topPane.add(targetDim);
		return topPane;
	}

	/*
	 * A method for creating the bottom panel and the submit button.
	 */
	private Component createBottom() {
		submitChoice = createButton("Submit", 's', false);
		submitChoice.setEnabled(false);
		bottomPane = new JPanel(new GridLayout(0, 1));
		bottomPane.add(submitChoice);
		return bottomPane;
	}

	/*
	 * A method that creates the right panel that contains the units based on
	 * the selected filters.
	 */
	private Component createRight() {
		for (int i = 0; i < allUnitStrings.length; i++) {
			selectedUnitsSet.addElement(allUnitStrings[i]);
		}
		allUnits = new JList(selectedUnitsSet);

		allUnits.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		allUnits.setVisibleRowCount(rowCount);
		allUnits.addListSelectionListener(new AllSelector());
		rightListPane = new JScrollPane(allUnits);
		return rightListPane;
	}

	/*
	 * A method to create the left panel that contains the used units.
	 */
	private Component createLeft() {
		usedUnits = new JList(usedUnitsModel);
		usedUnits.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		usedUnits.setVisibleRowCount(rowCount);
		usedUnits.addListSelectionListener(new UsedSelector());
		leftListPane = new JScrollPane(usedUnits);
		return leftListPane;
	}

	/*
	 * A method to create a button. Requires the displayed name of the button,
	 * the key for alt-key keyboard navigation, and a boolean value of its
	 * focusability.
	 */
	private JButton createButton(String title, char mnemonic, boolean focusable) {
		JButton currentButton = new JButton(title);
		currentButton.setMnemonic(mnemonic);
		currentButton.setFocusable(focusable);
		currentButton.addActionListener(acheck);
		return currentButton;
	}

	/*
	 * Another way to create a button that also includes an input of a boolean
	 * indicating if it starts off enabled.
	 */
	private JButton createButton(String title, char mnemonic,
			boolean focusable, boolean enabled) {
		JButton currentButton = createButton(title, mnemonic, focusable);
		currentButton.setEnabled(enabled);
		return currentButton;
	}

	/*
	 * A method to create a combo box from a string array, adding "All" to the
	 * top of the list.
	 */
	private JComboBox addComboBox(String[] stringArray) {
		String[] newStringArray = new String[stringArray.length + 1];
		newStringArray[0] = "All";
		for (int n = 0; n < stringArray.length; n++) {
			newStringArray[n + 1] = stringArray[n];
		}
		JComboBox currentBox = new JComboBox(newStringArray);
		currentBox.addActionListener(acheck);
		return currentBox;
	}

	/*
	 * A method that creates the center panel. It contains buttons for adding
	 * and removing units and changing the exponent each is raised to. It also
	 * contains the combo boxes that have the units systems and dimensions that
	 * can be used as filters.
	 */
	private Component createCenter() {

		// Create the buttons.
		addUnitButton = createButton("<", 'b', false);
		raisePower = createButton("^", 'p', false);
		removeUnitButton = createButton(">", 'f', false, false);
		lowerPower = createButton("v", 'n', false);

		// Create the combo boxes.
		dimensionChoice = addComboBox(allDimensionStrings);
		unitSystemChoice = addComboBox(allUnitSystemStrings);

		// Create the center panel.
		centerPane = new JPanel(new GridLayout(0, 1));

		// Add buttons and combo boxes to the center panel.
		centerPane.add(addUnitButton);
		centerPane.add(removeUnitButton);
		centerPane.add(raisePower);
		centerPane.add(lowerPower);
		centerPane.add(dimensionChoice);
		centerPane.add(unitSystemChoice);

		// Add borders to the panel.
		centerPane.setBorder(BorderFactory.createEmptyBorder(30, // top
				50, // left
				30, // bottom
				50) // right
				);

		return centerPane;
	}

	/*
	 * This creates a JFrame, adds the top, bottom, left, center, and right
	 * panels to it, and makes it appear.
	 */
	public void createAndShowGUI() {
		// JFrame.setDefaultLookAndFeelDecorated(true);
		JFrame frame = new JFrame("Unit Selector");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		UnitGraphics app = new UnitGraphics();
		app.setDim(targetDimension);
		Component centerContents = app.createCenter();
		Component leftContents = app.createLeft();
		Component rightContents = app.createRight();
		Component topContents = app.createTop();
		Component bottomContents = app.createBottom();
		frame.getContentPane().add(centerContents, BorderLayout.CENTER);
		frame.getContentPane().add(leftContents,
				BorderLayout.BEFORE_LINE_BEGINS);
		frame.getContentPane().add(rightContents, BorderLayout.AFTER_LINE_ENDS);
		frame.getContentPane().add(topContents, BorderLayout.BEFORE_FIRST_LINE);
		frame.getContentPane()
				.add(bottomContents, BorderLayout.AFTER_LAST_LINE);
		// Display the window.

		frame.pack();
		frame.setVisible(true);

	}

	/*
	 * A class that listens for changes to the selected unit.
	 */
	private class AllSelector implements ListSelectionListener {
		public void valueChanged(ListSelectionEvent event) {
			if (!event.getValueIsAdjusting() && !allUnits.isSelectionEmpty())

				selectedUnit = allUnits.getSelectedValue().toString();
			if (unitsVector.contains(selectedUnit)) {
				addUnitButton.setEnabled(false);
			} else {
				addUnitButton.setEnabled(true);
			}

		}
	}

	/*
	 * A class that listens for changes to the selected unit on the used unit
	 * list.
	 */
	private class UsedSelector implements ListSelectionListener {
		public void valueChanged(ListSelectionEvent event) {
			if (!event.getValueIsAdjusting() && !usedUnits.isSelectionEmpty())

				selectedUsedUnit = unitsVector.elementAt(
						usedUnits.getSelectedIndex()).toString();

		}
	}

	/*
	 * A class that listens for the non-list actions, including the button
	 * presses and the combo box selections.
	 */
	private class ActionCheck implements ActionListener {

		public void actionPerformed(ActionEvent e) {

			/*
			 * Create a string containing the current action command and display
			 * it in the top panel's label.
			 */
			String actionCommand = e.getActionCommand();
			label.setText("Action: " + actionCommand);

			if (actionCommand.equals("comboBoxChanged") && (switching == false)) {
				comboWasChanged(e);
			}

			if (actionCommand.equals("<") || actionCommand.equals(">")) {
				unitButtonPressed(e);
				checkUnitVsTarget(targetDimension);
			}
			if (actionCommand.equals("^") || actionCommand.equals("v")) {
				exponentButtonPressed(e);
				checkUnitVsTarget(targetDimension);
			}

			if (actionCommand.equals("Submit")) {
				int usize = unitsVector.size();
				Unit[] allUnits = new Unit[usize];
				double[] allExponents = new double[usize];
				for (int i = 0; i < usize; i++) {
					allUnits[i] = UnitFilter.stringToUnit(unitsVector
							.elementAt(i).toString());
					allExponents[i] = Double.parseDouble((exponentVector
							.elementAt(i).toString()));
				}
				currentCompoundUnit = new CompoundUnit(allUnits, allExponents);
				System.exit(0);
			}

		}

		private void checkUnitVsTarget(Dimension target) {
			int usize = unitsVector.size();
			Unit[] allUnits = new Unit[usize];
			double[] allExponents = new double[usize];
			for (int i = 0; i < usize; i++) {
				allUnits[i] = UnitFilter.stringToUnit(unitsVector.elementAt(i)
						.toString());
				allExponents[i] = Double.parseDouble(exponentVector
						.elementAt(i).toString());
			}
			currentCompoundUnit = new CompoundUnit(allUnits, allExponents);

			// System.out.println(currentCompoundUnit.dimension());

			if (currentCompoundUnit.dimension().equals(target)) {
				submitChoice.setEnabled(true);
			} else {
				submitChoice.setEnabled(false);
			}
		}

		private void exponentButtonPressed(ActionEvent e) {
			String actionCommand = e.getActionCommand();

			/*
			 * This is performed if the raise power button is pressed, there is
			 * a selected used unit.
			 */
			if (actionCommand.equals("^") && !selectedUsedUnit.equals(null)) {
				// The exponent and unit are stored.
				dl = exponentVector.elementAt(
						unitsVector.indexOf(selectedUsedUnit)).toString();
				di = Integer.valueOf(dl).intValue();
				// The exponent is increased by one.
				di++;
				// Some code to avoid having exponents be zero.
				if (di == 0) {
					di++;
				}
				// The values of the exponent and the display are updated.
				exponentVector.set(unitsVector.indexOf(selectedUsedUnit),
						String.valueOf(di));
				usedUnitsModel.set(unitsVector.indexOf(selectedUsedUnit),
						selectedUsedUnit + "^" + di);
			}

			/*
			 * This is performed if the lower power button is pressed.
			 */
			if (actionCommand.equals("v") && !selectedUsedUnit.equals(null)) {
				// The exponent and unit are stored.
				dl = exponentVector.elementAt(
						unitsVector.indexOf(selectedUsedUnit)).toString();
				di = Integer.valueOf(dl).intValue();
				// The exponent is lowered by one.
				di--;
				// Prevent exponent of zero by dropping it to -1.
				if (di == 0) {
					di--;
				}
				// The values of the exponent and the display are updated.
				exponentVector.set(unitsVector.indexOf(selectedUsedUnit),
						String.valueOf(di));
				usedUnitsModel.set(unitsVector.indexOf(selectedUsedUnit),
						selectedUsedUnit + "^" + di);
			}
		}

		private void unitButtonPressed(ActionEvent e) {
			String actionCommand = e.getActionCommand();

			/*
			 * This is excecuted if the command was to add the unit, there is a
			 * selected unit, and that unit is not yet used.
			 */
			if (actionCommand.equals("<")
					&& !unitsVector.contains(selectedUnit)
					&& !selectedUnit.equals(null)) {
				/*
				 * The unit is added to the units vector, and the default
				 * exponent of 1 is added to the exponent vector.
				 */
				unitsVector.addElement(selectedUnit);
				exponentVector.addElement("1");
				/*
				 * These two are combined with a "^" added to the dislayed used
				 * units list.
				 */
				usedUnitsModel.addElement(selectedUnit + "^"
						+ exponentVector.lastElement());
				// The remove unit button is enabled.
				removeUnitButton.setEnabled(true);
				// The added unit is selected on the used units list.
				usedUnits.setSelectedIndex(unitsVector.indexOf(selectedUnit));
				// The add unit button is disabled because the unit now in use.
				addUnitButton.setEnabled(false);
			}

			/*
			 * This is performed if the remove unit button hit and there is a
			 * selected used unit.
			 */
			if (actionCommand.equals(">") && !selectedUsedUnit.equals(null)
					&& !unitsVector.equals(null)) {

				/*
				 * Here is some code to decide which unit is selected after the
				 * current one is removed.
				 */
				di = unitsVector.indexOf(selectedUsedUnit);
				if (di > 0) {
					usedUnits.setSelectedIndex(di - 1);
					selectedUsedUnit = unitsVector.elementAt(di - 1).toString();

				} else {
					if (unitsVector.size() == 1) {
						removeUnitButton.setEnabled(false);
						usedUnits.clearSelection();
						selectedUsedUnit = null;
					} else {
						usedUnits.setSelectedIndex(di + 1);
						selectedUsedUnit = unitsVector.elementAt(di + 1)
								.toString();
					}
				}
				// Unit is removed from the display, unit, and exponent lists.
				usedUnitsModel.removeElementAt(di);
				unitsVector.removeElementAt(di);
				exponentVector.removeElementAt(di);
				/*
				 * The add unit button is enabled if there are units to add add
				 * the selected unit is not used.
				 */
				if (!selectedUnitsSet.equals(null)
						&& !unitsVector.contains(selectedUnit)) {
					addUnitButton.setEnabled(true);
				}
			}

		}

		private void comboWasChanged(ActionEvent e) {

			// If the unit system was changed...
			if (e.getSource().equals(unitSystemChoice)) {
				/*
				 * label.setText(unitSystemChoice.getSelectedItem()
				 * .toString());
				 */
				switching = true;
				dimensionChoice.setSelectedIndex(0);
				switching = false;
				selectedUnitsSet.removeAllElements();
				if (unitSystemChoice.getSelectedItem().equals("All")) {
					for (int i = 0; i < allUnitStrings.length; i++) {
						// System.out.println(newUnitSet[i]);
						selectedUnitsSet.addElement(allUnitStrings[i]);
					}
				} else {
					String[] unitList = Lister.uisString(unitSystemChoice
							.getSelectedItem().toString());
					for (int i = 0; i < unitList.length; i++) {
						//System.out.println(unitList[i]);
						selectedUnitsSet.addElement(unitList[i]);
					}
				}

				// If the dimension combo box was changed...
			} else if (e.getSource().equals(dimensionChoice)) {
				switching = true;
				unitSystemChoice.setSelectedIndex(0);
				switching = false;
				if (dimensionChoice.getSelectedItem().toString() == "All") {
					selectedUnitsSet.removeAllElements();
					for (int i = 0; i < allUnitStrings.length; i++) {
						// System.out.println(newUnitSet[i]);
						selectedUnitsSet.addElement(allUnitStrings[i]);
					}
				} else {
					label.setText(dimensionChoice.getSelectedItem().toString());
					String[] newUnitSet = UnitFilter.filter(UnitFilter
							.stringToDim(dimensionChoice.getSelectedItem()
									.toString()), allUnitStrings);
					selectedUnitsSet.removeAllElements();
					for (int i = 0; i < newUnitSet.length; i++) {
						// System.out.println(newUnitSet[i]);
						selectedUnitsSet.addElement(newUnitSet[i]);
					}
					if (selectedUnitsSet.equals(null)) {
						addUnitButton.setEnabled(false);
					} else {
						addUnitButton.setEnabled(true);

					}

				}
			}
		}

	}

}
