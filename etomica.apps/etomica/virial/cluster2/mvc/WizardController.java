package etomica.virial.cluster2.mvc;

/**
 * Simple wizard template.
 *
 * @author Demian Lessa
 */
public abstract class WizardController implements Runnable, ViewResponseListener {

  private WizardView wizard;
  private State internalState = new DefaultState();

  /**
   * Invoked later by the Swing thread.
   *
   */
  public void run() {

    // initialize the state, create the wizard
    preProcess();
    // first wizard page
    WizardPageView pageView = nextPageView(null, null);
    // CONTRACT: there must exist a view corresponding to the action response
    assert (pageView != null);
    // the wizard attaches its public UI to the page view
    getWizard().attachPageView(pageView);
    // configure the page view with the internal state and the listener
    pageView.configure(getState());
    pageView.setResponseListener(this);
    pageView.display();
    // display the wizard; after this first display, all other displays are event-based
    // and are indirectly controlled by the event handlers in each page view
    getWizard().display();
  }

  // Returns the actual wizard on top of which the wizard pages will be displayed
  protected WizardView getWizard() {

    if (wizard == null) {
      createWizard();
    }
    return wizard;
  }

  public void onViewResponse(ViewResponse viewResponse) {

    // CONTRACT: the wizard page must return a valid response
    assert (viewResponse != null);
    Action action = nextAction(viewResponse);
    // CONTRACT: there must exist an action corresponding to the view response
    assert (action != null);
    ActionResponse actionResponse = action.execute(viewResponse, getState());
    // CONTRACT: every action must return a response object
    assert (actionResponse != null);
    // detach the current view because we are done with it
    getWizard().detachPageView((WizardPageView) viewResponse.getView());
    if (!actionResponse.getStatus().isTerminated()) {
      WizardPageView pageView = nextPageView(actionResponse, viewResponse);
      // CONTRACT: there must exist a view corresponding to the action response
      assert (pageView != null);
      // the wizard attaches its public UI to the page view
      getWizard().attachPageView(pageView);
      pageView.configure(getState());
      pageView.setResponseListener(this);
      pageView.display();
    }
    else {
      // post-processing with the last action response
      postProcess(actionResponse);
      getWizard().close();
    }
  }

  /**
   * Performs all necessary initialization prior to the wizard execution. A convenience
   * default implementation is provided to initialize the state and create the wizard.
   */
  protected void preProcess() {

    wizard = createWizard();
    initializeState();
    getWizard().configure(getState());
  }

  /**
   * Initializes the wizard view container.
   */
  protected abstract WizardView createWizard();

  /**
   * Initializes the internal state to be used by the wizard. The default
   * implementation is a no-operation.
   */
  protected void initializeState() {

    // no-operation
  }

  /**
   * Returns this workflow's state instance.
   */
  protected State getState() {

    return internalState;
  }

  /**
   * Performs any necessary result processing after the wizard has completed. The default
   * implementation is a no-operation.
   *
   */
  protected void postProcess(ActionResponse actionResponse) {

    // no-operation
  }

  /**
   * Find and instantiate the next action to execute based on the current state of the
   * wizard and the last view response.
   */
  protected abstract Action nextAction(ViewResponse response);

  /**
   * Find and instantiate the next view to display based on the current state of he wizard
   * and the last action response.
   *
   */
  protected abstract WizardPageView nextPageView(ActionResponse actionResponse, ViewResponse viewResponse);
}