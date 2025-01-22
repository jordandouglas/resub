package resub.substitutionmodel;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.core.Input.Validate;
import beast.base.core.Log;

@Description("Initialises indicators for the stochastic variable selector substitition model")
public class SVSInit extends BEASTObject implements StateNodeInitialiser {
	
	
	public Input<BooleanParameter> indicatorInput = new Input<BooleanParameter>("indicator", "boolean indicator parameter", Validate.REQUIRED);
	public Input<List<Alignment>> dataInput = new Input<List<Alignment>>("data", "alignment to use for initialising indicator", new ArrayList<>());
	public Input<DataType> dataTypeInput = new Input<>("dataType", "the data type", Input.Validate.REQUIRED);
	
	int nrOfStates;
	boolean isSymmetric;
	
	@Override
	public void initAndValidate() {
		
		
		if (dataInput.get().isEmpty()) {
			throw new IllegalArgumentException("Please provide at least 1 data");
		}
		
		
		// Check dimensions
		this.nrOfStates = dataTypeInput.get().getStateCount();
		if (indicatorInput.get().getDimension() == nrOfStates*(nrOfStates-1)) {
			this.isSymmetric = false;
		}
		
		else if (indicatorInput.get().getDimension() == nrOfStates*(nrOfStates-1)/2) {
			this.isSymmetric = true;
		}
		
		else {
			throw new IllegalArgumentException("inidicator should have either " + (nrOfStates*(nrOfStates-1)) + " or " + (nrOfStates*(nrOfStates-1)/2) + " dimensions but it has " + indicatorInput.get().getDimension());
		}
		
		
		
	}

	@Override
	public void initStateNodes() {
		
		// If there are no occurrences of si and sj coexisting in the same column, set its indicator to false
		BooleanParameter indicator = indicatorInput.get();
		boolean value = true;
		int count = 0;
		for (int i = 0; i < nrOfStates; i++) {
			
			// Symmetric
			if (this.isSymmetric) {
				
				
				for (int j = i+1; j <  nrOfStates; j++) {
					value = statesCoexist(i, j, dataInput.get(), dataTypeInput.get());
					indicator.setValue(count, value);
					count++;
				}
				
				
			}
			
			
			// Asymmetric
			else {
				
	            for (int j = 0; j < i; j++) {
	            	value = statesCoexist(i, j, dataInput.get(), dataTypeInput.get());
	            	indicator.setValue(i * (nrOfStates - 1) + j, value);
	            }
	            for (int j = i + 1; j < nrOfStates; j++) {
	            	value = statesCoexist(i, j, dataInput.get(), dataTypeInput.get());
	            	indicator.setValue(i * (nrOfStates - 1) + j - 1, value);
	            }
				
			}
			
			
		}
		
	
		
	}
	
	
	/**
	 * Does there exist a site where the two states both exist?
	 * @param state1
	 * @param state2
	 * @param data
	 * @return
	 */
	public static boolean statesCoexist(int state1, int state2, List<Alignment> alignments, DataType dt) {
		
		for (Alignment data : alignments) {
			
			for (int patternNum = 0; patternNum < data.getPatternCount(); patternNum++) {
				
				
				// Check if this pattern has both state1 and state2
				int[] codes = data.getPattern(patternNum);
				boolean has1 = false;
				boolean has2 = false;
				for (int taxonNr = 0; taxonNr < codes.length; taxonNr++) {
					int code = codes[taxonNr];
					if (code == state1) has1 = true;
					if (code == state2) has2 = true;
				}
				
				if (has1 && has2) {
					//Log.warning(dt.getCharacter(state1) + " and " + dt.getCharacter(state2) + " DO coexist in " + data.getID() + " as " + state1 + " and " + state2);
					return true;
				}
				
			}
			
			//Log.warning(dt.getCharacter(state1) + " and " + dt.getCharacter(state2) + " do not coexist in " + data.getID());
			
			
		}
		
	
		
		//Log.warning(dt.getCharacter(state1) + " and " + dt.getCharacter(state2) + " do not coexist");
		return false;
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(indicatorInput.get());
	}

	
	
}








