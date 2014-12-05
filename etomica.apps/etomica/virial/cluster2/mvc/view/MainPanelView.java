/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.BorderLayout;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JSplitPane;

import com.jgoodies.forms.factories.Borders;
import com.jgoodies.uif_lite.component.Factory;

public class MainPanelView {

  public static JComponent build() {

    JPanel panel = new JPanel(new BorderLayout());
    panel.setOpaque(false);
    panel.setBorder(Borders.DIALOG_BORDER);
    panel.add(buildHorizontalSplit());
    return panel;
  }

  private static JComponent buildHorizontalSplit() {

    JSplitPane pane = Factory.createStrippedSplitPane(
        JSplitPane.HORIZONTAL_SPLIT, LeftPanelView.build(), RightPanelView
            .build(), 0.2f);
    pane.setOpaque(false);
    return pane;
  }
}
