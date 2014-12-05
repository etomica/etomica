/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.AbstractButton;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.KeyStroke;

import com.jgoodies.looks.Options;
import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.windows.WindowsLookAndFeel;

public class ToolBarView {

  /*
   * Builds, configures and returns the Requests HeaderStyle, look-specific
   * BorderStyles, and Plastic 3D Hint from Launcher.
   */
  public static JToolBar mainToolbar(ApplicationView parent) {

    JGoodiesSettings settings = ApplicationUI.uiSettings;
    JToolBar toolBar = new JToolBar();
    toolBar.setFloatable(false);
    toolBar.putClientProperty("JToolBar.isRollover", Boolean.TRUE);
    toolBar.putClientProperty(Options.HEADER_STYLE_KEY, settings
        .getToolBarHeaderStyle());
    toolBar.putClientProperty(PlasticLookAndFeel.BORDER_STYLE_KEY, settings
        .getToolBarPlasticBorderStyle());
    toolBar.putClientProperty(WindowsLookAndFeel.BORDER_STYLE_KEY, settings
        .getToolBarWindowsBorderStyle());
    toolBar.putClientProperty(PlasticLookAndFeel.IS_3D_KEY, settings
        .getToolBar3DHint());
    AbstractButton button;
    toolBar.add(createToolBarButton("backward.gif", "Back"));
    button = createToolBarButton("forward.gif", "Next");
    button.setEnabled(false);
    toolBar.add(button);
    toolBar.add(createToolBarButton("home.gif", "Home"));
    toolBar.addSeparator();
    ActionListener openAction = parent.new OpenFileActionListener();
    button = createToolBarButton("open.gif", "Open", openAction, KeyStroke
        .getKeyStroke(KeyEvent.VK_O, KeyEvent.CTRL_DOWN_MASK));
    button.addActionListener(openAction);
    toolBar.add(button);
    toolBar.add(createToolBarButton("print.gif", "Print"));
    toolBar.add(createToolBarButton("refresh.gif", "Update"));
    toolBar.addSeparator();
    ButtonGroup group = new ButtonGroup();
    button = createToolBarRadioButton("pie_mode.png", "Pie Chart");
    button
        .setSelectedIcon(ApplicationUI.readImageIcon("pie_mode_selected.gif"));
    group.add(button);
    button.setSelected(true);
    toolBar.add(button);
    button = createToolBarRadioButton("bar_mode.png", "Bar Chart");
    button
        .setSelectedIcon(ApplicationUI.readImageIcon("bar_mode_selected.gif"));
    group.add(button);
    toolBar.add(button);
    button = createToolBarRadioButton("table_mode.png", "Table");
    button.setSelectedIcon(ApplicationUI
        .readImageIcon("table_mode_selected.gif"));
    group.add(button);
    toolBar.add(button);
    toolBar.addSeparator();
    button = createToolBarButton("help.gif", "Open Help");
    button.addActionListener(parent.new HelpActionListener());
    toolBar.add(button);
    return toolBar;
  }

  /**
   * Creates and returns a JButton configured for use in a J
   * <p>
   * This is a simplified method that is overridden by the Looks Demo. The full
   * code uses the JGoodies UI framework's ToolBarButton that better handles
   * platform differences.
   */
  private static AbstractButton createToolBarButton(String iconName,
      String toolTipText) {

    JButton button = new JButton(ApplicationUI.readImageIcon(iconName));
    button.setToolTipText(toolTipText);
    button.setFocusable(false);
    return button;
  }

  private static AbstractButton createToolBarButton(String iconName,
      String toolTipText, ActionListener action, KeyStroke keyStroke) {

    AbstractButton button = createToolBarButton(iconName, toolTipText);
    button.registerKeyboardAction(action, keyStroke,
        JComponent.WHEN_IN_FOCUSED_WINDOW);
    return button;
  }

  /**
   * Creates and returns a JToggleButton configured for use in a J
   * <p>
   * This is a simplified method that is overridden by the Looks Demo. The full
   * code uses the JGoodies UI framework's ToolBarButton that better handles
   * platform differences.
   */
  private static AbstractButton createToolBarRadioButton(String iconName,
      String toolTipText) {

    JToggleButton button = new JToggleButton(ApplicationUI
        .readImageIcon(iconName));
    button.setToolTipText(toolTipText);
    button.setFocusable(false);
    return button;
  }
}