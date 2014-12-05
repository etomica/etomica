/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.event.*;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;

import etomica.virial.cluster2.mvc.*;

public class DefaultWizardPage implements WizardPageView {

  private ViewResponseListener viewResponseListener;
  private MouseAdapter mouseAdapter;
  private KeyAdapter keyAdapter;
  private WizardController controller;

  /**
   * Creates the wizard page and default keyboard and mouse handlers to delegate the
   * handling to the appropriate methods.
   */
  public DefaultWizardPage(WizardController controller) {

    this.controller = controller;
    keyAdapter = new KeyAdapter() {

      @Override
      public void keyPressed(KeyEvent e) {

        if (e.getComponent().isEnabled()) {
          if (e.getKeyCode() == KeyEvent.VK_SPACE || e.getKeyCode() == KeyEvent.VK_ENTER) {
            if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_NEXT_BUTTON)) {
              doNext();
            }
            else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_BACK_BUTTON)) {
              doBack();
            }
            else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_CANCEL_BUTTON)) {
              doCancel();
            }
            else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_HELP_BUTTON)) {
              doHelp();
            }
            else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_FINISH_BUTTON)) {
              doFinish();
            }
          }
        }
      }
    };

    mouseAdapter = new MouseAdapter() {

      @Override
      public void mouseClicked(MouseEvent e) {

        if (e.getComponent().isEnabled()) {
          if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_NEXT_BUTTON)) {
            doNext();
          }
          else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_BACK_BUTTON)) {
            doBack();
          }
          else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_CANCEL_BUTTON)) {
            doCancel();
          }
          else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_HELP_BUTTON)) {
            doHelp();
          }
          else if (e.getComponent() == getController().getState().getProperty(ClusterWizard.KEY_FINISH_BUTTON)) {
            doFinish();
          }
        }
      }
    };
  }

  public WizardController getController() {

    return controller;
  }

  public void attach(String key, Object object) {

    getController().getState().setProperty(key, object);

    if (object instanceof JButton) {
      attachButton(key, (JButton) object);
    }
    else if (object instanceof JPanel) {
      attachPanel(key, (JPanel) object);
    }
  }

  protected void attachButton(String key, JButton button) {

    button.addMouseListener(mouseAdapter);
    button.addKeyListener(keyAdapter);
  }

  public void attachDone() {

    loadFromState();
  }

  protected void attachPanel(String key, JPanel panel) {

    // default no-op attach
  }

  protected ViewResponse createResponse(final ViewStatus status, final List<MVCException> errors) {

    final WizardPageView thisView = this;
    return new ViewResponse() {

      public List<MVCException> getErrors() {

        return errors;
      }

      public ViewStatus getStatus() {

        return status;
      }

      public View getView() {

        return thisView;
      }

    };
  }

  public void detach(String key, Object object) {

    if (object instanceof JButton) {
      detachButton(key, (JButton) object);
    }
    if (key.startsWith(ClusterWizard.KEY_PANE)) {
      detachPanel(key, (JPanel) object);
    }
  }

  protected void detachButton(String key, JButton button) {

    button.removeMouseListener(mouseAdapter);
    button.removeKeyListener(keyAdapter);
  }

  public void detachDone() {

    display();
  }

  protected void detachPanel(String key, JPanel panel) {

    panel.removeAll();
  }

  public void display() {

    JPanel modelPane = (JPanel) getController().getState().getProperty(ClusterWizard.KEY_MODEL_PANE);
    modelPane.validate();
    modelPane.repaint();
    JPanel figurePane = (JPanel) getController().getState().getProperty(ClusterWizard.KEY_FIGURE_PANE);
    figurePane.validate();
    figurePane.repaint();
  }

  protected void doBack() {

    // this is conservative but a good idea since state on a page X may depend on the state of a page X-1
    rollbackChanges();
    viewResponseListener.onViewResponse(createResponse(ViewStatus.CONTINUE_PRIOR, null));
  }

  protected void doCancel() {

    commitChanges();
    viewResponseListener.onViewResponse(createResponse(ViewStatus.COMPLETE_USER_CANCELED, null));
  }

  protected void doFinish() {

    commitChanges();
    viewResponseListener.onViewResponse(createResponse(ViewStatus.COMPLETE_SUCCESS, null));
  }

  protected void doHelp() {

    viewResponseListener.onViewResponse(createResponse(ViewStatus.CONTINUE_HELP, null));
  }

  protected void doNext() {

    commitChanges();
    viewResponseListener.onViewResponse(createResponse(ViewStatus.CONTINUE_NEXT, null));
  }

  public void setResponseListener(ViewResponseListener listener) {

    this.viewResponseListener = listener;
  }

  public void loadFromState() {

    // default no-op load
  }

  public void commitChanges() {

    // default no-op commit
  }

  public void rollbackChanges() {

    // sets the state of this wizard page back to its default
    getController().getState().loadDefaultState(getPageId());
  }

  public void initializeUI() {

    // default no-op initialize
  }

  public int getPageId() {

    return -1;
  }
}