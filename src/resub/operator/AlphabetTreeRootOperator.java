package resub.operator;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONStringer;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.StateNode;
import resub.substitutionmodel.epochs.AlphabetEpochs;
import beast.base.evolution.tree.Tree;

/**
 * 
 * @author Jordan Douglas
 */ 
@Description("An operator that calls another operator but then adjusts its hastings ratio if the tree root has changed height to reflect changes in the alphabet epoch ages")
public class AlphabetTreeRootOperator extends Operator {
	

    final public Input<Tree> treeInput = new Input<>("tree", "the tree that is being changed", Input.Validate.REQUIRED);
    final public Input<Operator> operatorInput = new Input<>("operator", "the tree operator", Input.Validate.REQUIRED);
    final public Input<AlphabetEpochs> alphabetInput = new Input<>("alphabet", "the alphabet ofepochs, so that we know how many alphabets there are", Input.Validate.REQUIRED);

	@Override
	public void initAndValidate() {

		

		
	}

	@Override
	public double proposal() {
		
		Operator operator = operatorInput.get();
		double rootHeightBefore = treeInput.get().getRoot().getHeight();
		double logHR = operator.proposal();
		double rootHeightAfter = treeInput.get().getRoot().getHeight();
		
		// Adjust Jacobian
		if (rootHeightAfter != rootHeightBefore) {
			int nalphabets = alphabetInput.get().getNOldEpochs();
			logHR += nalphabets * (Math.log(rootHeightAfter) - Math.log(rootHeightBefore));
			//Log.warning(rootHeightBefore + " to " + rootHeightAfter);
		}
		
		return logHR;

	}
	
	
	
	 @Override
	    public void setOperatorSchedule(final OperatorSchedule operatorSchedule) {
	    	super.setOperatorSchedule(operatorSchedule);
	    	Operator operator = operatorInput.get();
	    	operator.setOperatorSchedule(operatorSchedule);
	    }
	
	@Override
	public void accept() {
		Operator operator = operatorInput.get();
		operator.accept();
		super.accept();
	}
	
	
	@Override
	public void reject(int reason) {
		Operator operator = operatorInput.get();
		operator.reject(reason);
		super.reject(reason);
	}
	
	
    @Override
    public void optimize(double logAlpha) {
    	Operator operator = operatorInput.get();
    	operator.optimize(logAlpha);
    }
	

    
    @Override
    public List<StateNode> listStateNodes() {
    	Operator operator = operatorInput.get();
    	return operator.listStateNodes();
    }
    

    
    
    
    @Override
    public void storeToFile(final PrintWriter out) {
    	

        
    	try {
	        JSONStringer json = new JSONStringer();
	        json.object();
	
	        if (getID() == null) setID("unknown");
	
	        // id
	        json.key("id").value(getID());
	        
	        Operator operator = operatorInput.get();
	        
	        
	        
	        // Store generic beast core properties by writing its json to a string and then parsing it back
	        StringWriter outStr = new StringWriter();
	        PrintWriter writer = new PrintWriter(outStr);
	        super.storeToFile(writer);
	        JSONObject obj = new JSONObject(outStr.toString());
	        for (String key : obj.keySet()) {
	        	if (key.equals("id")) continue;
	        	json.key(key).value(obj.get(key));
	        }
	        
	
	        // Store sub-operator in a list
	        JSONArray operatorListJson = new JSONArray();
	        outStr = new StringWriter();
	        writer = new PrintWriter(outStr);
        	operator.storeToFile(writer);
        	obj = new JSONObject(outStr.toString());
        	operatorListJson.put(obj);
	        json.key("operator").value(operatorListJson);
	        
	        json.endObject();
	        out.print(json.toString());
	        

	        
    	} catch (JSONException e) {
    		// failed to log operator in state file
    		// report and continue
    		e.printStackTrace();
    	}
    }
    
    
    

    @Override
    public void restoreFromFile(JSONObject o) {

    	
    	super.restoreFromFile(o);
    	
    	
    	try {
    		
	    	// Load sub-operator
	        JSONArray operatorArray = o.getJSONArray("operator");
	        if (operatorArray.length() != 1) {
	        	return;
        	}
	        
	        Operator operator = operatorInput.get();
        	JSONObject jsonOp = operatorArray.getJSONObject(0);
	        operator.restoreFromFile(jsonOp);
    		
	        super.restoreFromFile(o);  	
	        
    	} catch (JSONException e) {
    		// failed to restore from state file
    		// report and continue
    		e.printStackTrace();
    	}
    }

    
    
    

}







