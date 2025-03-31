package resub.operator;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.util.Randomizer;
import resub.substitutionmodel.epochs.AlphabetEpoch;
import resub.substitutionmodel.epochs.AlphabetEpochs;


@Description("Samples the alpha and beta states of an epoch, in a way that guarantees the epoch does not offend any earlier epochs")
public class EpochStateUniform extends Operator {
	
	final public Input<AlphabetEpochs> epochsInput = new Input<>("epochs", "vector of alphabet epochs", Input.Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		
		
	}

	@Override
	public double proposal() {
		
		AlphabetEpochs epochs = epochsInput.get();
		
		// Sample an epoch
		epochs.update();
		int epochIndex = Randomizer.nextInt(epochs.getNOldEpochs());
		AlphabetEpoch epoch = epochs.getEpochByIndex(epochIndex);
		int nstates = epochsInput.get().getNStatesTotal();
		
		if (nstates < 2) {
			Log.warning("Unexpected: there are " + nstates + " in epoch " + epochIndex);
			return Double.NEGATIVE_INFINITY;
		}
		
		
		// Sample new alpha and beta uniformly
		int alphaNew = Randomizer.nextInt(nstates);
		int betaNew = Randomizer.nextInt(nstates);
		
		
		// Check that alpha and beta are not pointing to invalid states
		while (!isValid(alphaNew, betaNew, nstates, epoch)) {
			alphaNew = Randomizer.nextInt(nstates);
			betaNew = Randomizer.nextInt(nstates);
		}
		
		
		epochs.getAlphas().setValue(epochIndex, alphaNew);
		epochs.getBetas().setValue(epochIndex, betaNew);
		

		
		//Log.warning(nstates + " -> " + alphaNew + " " + betaNew);
		
		return 0;
	}
	
	@Override
	public List<StateNode> listStateNodes() {
		List<StateNode> nodes = new ArrayList<>();
		nodes.add(epochsInput.get().getAlphas());
		nodes.add(epochsInput.get().getBetas());
		return nodes;
	}
	
	
	private boolean isValid(int alpha, int beta, int nstates, AlphabetEpoch epoch) {
		
		if (alpha == beta) return false;
		
		String stateNameInEpoch = epoch.getStateNames().get(alpha);
		if (stateNameInEpoch.equals(AlphabetEpoch.BLANK_STATE)) {
			return false;
		}
		
		
		return true;
		
		
	}

}
