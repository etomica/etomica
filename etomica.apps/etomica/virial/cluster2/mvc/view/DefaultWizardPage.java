package etomica.virial.cluster2.mvc.view;

import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.KeyAdapter;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;

import etomica.virial.cluster2.mvc.DefaultState;
import etomica.virial.cluster2.mvc.MVCException;
import etomica.virial.cluster2.mvc.State;
import etomica.virial.cluster2.mvc.View;
import etomica.virial.cluster2.mvc.ViewResponse;
import etomica.virial.cluster2.mvc.ViewResponseListener;
import etomica.virial.cluster2.mvc.ViewStatus;
import etomica.virial.cluster2.mvc.WizardPageView;

public class DefaultWizardPage implements WizardPageView {

  private ViewResponseListener viewResponseListener;
  private State data = new DefaultState();
  private MouseAdapter mouseAdapter;
  private KeyAdapter keyAdapter;

  /**
   * Creates the wizard page and default keyboard and mouse handlers to delegate the
   * handling to the appropriate methods.
   */
  public DefaultWizardPage() {

    keyAdapter = new KeyAdapter() {

      @Override
      public void keyPressed(KeyEvent e) {

        if (e.getComponent().isEnabled()) {
          if (e.getKeyCode() == KeyEvent.VK_SPACE || e.getKeyCode() == KeyEvent.VK_ENTER) {
            if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_NEXT_BUTTON)) {
              doNext();
            }
            else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_BACK_BUTTON)) {
              doBack();
            }
            else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_CANCEL_BUTTON)) {
              doCancel();
            }
            else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)) {
              doHelp();
            }
            else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_FINISH_BUTTON)) {
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
          if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_NEXT_BUTTON)) {
            doNext();
          }
          else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_BACK_BUTTON)) {
            doBack();
          }
          else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_CANCEL_BUTTON)) {
            doCancel();
          }
          else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)) {
            doHelp();
          }
          else if (e.getComponent() == getData().getProperty(ClusterWizard.KEY_FINISH_BUTTON)) {
            doFinish();
          }
        }
      }
    };
  }

  public void attach(String key, Object object) {

    getData().setProperty(key, object);

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

    // no-op here; descendant classes can override this
  }

  protected void attachPanel(String key, JPanel panel) {

    // no-op here; descendant classes can override this
  }

  public void configure(State state) {

    for (String key : state.getKeys()) {
      getData().setProperty(key, state.getProperty(key));
    }
  }

  protected ViewResponse createResponse(final ViewStatus status, final List<MVCException> errors) {

    final WizardPageView thisView = this;
    return new ViewResponse() {

      public State getData() {

        return getData();
      }

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

    JPanel modelPane = (JPanel) getData().getProperty(ClusterWizard.KEY_MODEL_PANE);
    modelPane.validate();
    modelPane.repaint();
    JPanel figurePane = (JPanel) getData().getProperty(ClusterWizard.KEY_FIGURE_PANE);
    figurePane.validate();
    figurePane.repaint();
  }

  protected void doBack() {

    viewResponseListener.onViewResponse(createResponse(ViewStatus.CONTINUE_PRIOR, null));
  }

  protected void doCancel() {

    viewResponseListener.onViewResponse(createResponse(ViewStatus.COMPLETE_USER_CANCELED, null));
  }

  protected void doFinish() {

    viewResponseListener.onViewResponse(createResponse(ViewStatus.COMPLETE_SUCCESS, null));
  }

  protected void doHelp() {

    viewResponseListener.onViewResponse(createResponse(ViewStatus.CONTINUE_HELP, null));
  }

  protected void doNext() {

    viewResponseListener.onViewResponse(createResponse(ViewStatus.CONTINUE_NEXT, null));
  }

  protected State getData() {

    return data;
  }

  public void setResponseListener(ViewResponseListener listener) {

    this.viewResponseListener = listener;
  }
}