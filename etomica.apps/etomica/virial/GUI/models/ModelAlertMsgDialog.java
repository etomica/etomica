package etomica.virial.GUI.models;

public class ModelAlertMsgDialog {
	
	private String alertMessage;
	
	private ModelSpeciesSelection modelSpeciesSelection;
	
	public ModelAlertMsgDialog(String alertMessage){
		setAlertMessage(alertMessage);
	}

	public ModelSpeciesSelection getSpeciesSelectionDM() {
		return modelSpeciesSelection;
	}

	public void setSpeciesSelectionDM(ModelSpeciesSelection modelSpeciesSelection) {
		this.modelSpeciesSelection = modelSpeciesSelection;
	}

	public String getAlertMessage() {
		return this.alertMessage;
	}

	public void setAlertMessage(String alertMessage) {
		this.alertMessage = alertMessage;
	}
	
}
