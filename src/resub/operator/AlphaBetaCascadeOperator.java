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


@Description("Samples the alpha and beta states from epoch i through to the oldest epoch, in a way that guarantees the resulting state is valid")
public class AlphaBetaCascadeOperator extends Operator {
	
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
		
		int nstates = epochsInput.get().getNStatesTotal();
		if (nstates < 2) {
			Log.warning("Unexpected: there are " + nstates + " states");
			return Double.NEGATIVE_INFINITY;
		}
		
		
		for (int epochNr = epochIndex; epochNr < epochs.getNOldEpochs(); epochNr++) {
			
			epochs.requestUpdate();
			epochs.update();
			AlphabetEpoch epoch = epochs.getEpochByRank(epochNr);
			
			// Sample new alpha and beta uniformly
			int alphaNew = Randomizer.nextInt(nstates);
			int betaNew = Randomizer.nextInt(nstates);
			
			
			// Check that alpha and beta are not pointing to invalid states
			while (!isValid(alphaNew, betaNew, nstates, epoch)) {
				alphaNew = Randomizer.nextInt(nstates);
				betaNew = Randomizer.nextInt(nstates);
			}
			
			
			epochs.getAlphas().setValue(epochNr, alphaNew);
			epochs.getBetas().setValue(epochNr, betaNew);
			
			
			
		}
		
		//Log.warning(nstates + " -> " + alphaNew + " " + betaNew);
		
		return 0; // Is this correct??
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
