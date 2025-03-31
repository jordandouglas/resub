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


@Description("Swaps the alpha and beta of an epoch")
public class AlphaBetaSwapper extends Operator {
	
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
		int alphaOld = epochs.getAlphas().getValue(epochIndex);
		int betaOld = epochs.getBetas().getValue(epochIndex);
		
		
		
		epochs.getAlphas().setValue(epochIndex, betaOld);
		epochs.getBetas().setValue(epochIndex, alphaOld);
		
		
		return 0;
	}
	
	@Override
	public List<StateNode> listStateNodes() {
		List<StateNode> nodes = new ArrayList<>();
		nodes.add(epochsInput.get().getAlphas());
		nodes.add(epochsInput.get().getBetas());
		return nodes;
	}


}
