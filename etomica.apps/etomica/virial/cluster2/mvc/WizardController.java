/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

/**
 * Simple wizard template.
 *
 * @author Demian Lessa
 */
public abstract class WizardController implements Runnable, ViewResponseListener {

  private WizardView view;
  private WizardState state;

  /**
   * Invoked later by the Swing thread.
   */
  public void run() {

    // initialize the state, create the wizard
    preProcess();
    // first wizard page
    WizardPageView pageView = createPageView(null, null);
    // prepare the state
    prepareState(pageView);
    // prepare the page view
    preparePageView(pageView);
    // display the page view
    pageView.display();
    // display the wizard; after this first display, all other displays are event-based
    // and are indirectly controlled by the event handlers in each page view
    getWizardView().display();
  }

  private void preparePageView(WizardPageView pageView) {

    // initialize the UI before actually displaying it
    pageView.initializeUI();
    // the wizard attaches its public UI to the page view
    getWizardView().attachPageView(pageView);
    // configure the page view with the listener
    pageView.setResponseListener(this);
  }

  private void prepareState(WizardPageView pageView) {

    // initialize the part of the state corresponding to this page, if necessary
    if (!getState().isStateLoaded(pageView.getPageId())) {
      getState().loadDefaultState(pageView.getPageId());
    }
  }

  // Returns the actual wizard on top of which the wizard pages are displayed
  protected WizardView getWizardView() {

    if (view == null) {
      createWizard();
    }
    return view;
  }

  public void onViewResponse(ViewResponse viewResponse) {

    // CONTRACT: the wizard page must return a valid response
    assert (viewResponse != null);
    Action action = nextAction(viewResponse);
    // CONTRACT: there must exist an action corresponding to the view response
    assert (action != null);
    ActionResponse actionResponse = action.execute(viewResponse);
    // CONTRACT: every action must return a response object
    assert (actionResponse != null);
    // detach the current view because we are done with it
    getWizardView().detachPageView((WizardPageView) viewResponse.getView());
    if (!actionResponse.getStatus().isTerminated()) {
      WizardPageView pageView = createPageView(actionResponse, viewResponse);
      // prepare the state
      prepareState(pageView);
      // prepare the page view
      preparePageView(pageView);
      // display the page view
      pageView.display();
    }
    else {
      // post-processing with the last action response
      postProcess(actionResponse);
      getWizardView().close();
    }
  }

  /**
   * Performs all necessary initialization prior to the wizard execution. A convenience
   * default implementation is provided to initialize the state and create the wizard.
   */
  protected void preProcess() {

    state = createWizardState();
    view = createWizard();
    getWizardView().initializeUI();
  }

  /**
   * Initializes the wizard view container.
   */
  protected abstract WizardView createWizard();

  /**
   * Create the actual wizard state instance.
   */
  protected abstract WizardState createWizardState();

  /**
   * Returns this workflow's state instance.
   */
  public WizardState getState() {

    return state;
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
  protected abstract WizardPageView createPageView(ActionResponse actionResponse, ViewResponse viewResponse);
}