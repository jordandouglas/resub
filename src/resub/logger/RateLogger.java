package resub.logger;

import java.io.PrintStream;

import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;

@Description("Logger for GTR rates that produces readable labels in trace log")
public class RateLogger extends BEASTObject implements Loggable {
	final public Input<RealParameter> parameterInput = new Input<>("parameter", "rate parameter representing rates of general substitution model", Validate.REQUIRED);
	final public Input<DataType> dataTypeInput = new Input<>("dataType", "data type for providing transition rate labels", Validate.OPTIONAL);

	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "boolean indicator for which rates are non-zero (optional)", Validate.OPTIONAL);
	
	
	
	private RealParameter parameter;
	private DataType datatype;
	
	private boolean symmetric;
	
	
	@Override
	public void initAndValidate() {
		parameter = parameterInput.get();
		datatype = dataTypeInput.get();
		
		if (datatype != null) {
		
			
			int expectedDimSymmetric = (datatype.getStateCount())*(datatype.getStateCount()-1) / 2;
			int expectedDimAsymmetric = (datatype.getStateCount())*(datatype.getStateCount()-1);
			
			if (parameter.getDimension() == expectedDimSymmetric) {
				this.symmetric = true;
			}
			
			else if (parameter.getDimension() == expectedDimAsymmetric) {
				this.symmetric = false;
			}
			
			else {
				throw new IllegalArgumentException("Expected either " + expectedDimSymmetric + " (symmetric) or " + expectedDimAsymmetric + " (asymmetric) rates but there are " + parameter.getDimension());
			}
			
		
		}
		
		
		if (indicatorInput.get() != null && indicatorInput.get().getDimension() != parameter.getDimension()) {
			throw new IllegalArgumentException("Expected " + parameter.getDimension() + " indicator values but there are " + indicatorInput.get().getDimension());
		}

		
	}

	@Override
	public void init(PrintStream out) {
	
		
		if (datatype != null) {
		
			for (int i = 0; i < datatype.getStateCount(); i++) {
			
				if (this.symmetric) {
					
					// Symmetric
					for (int j = i+1; j < datatype.getStateCount(); j++) {
						out.append(((BEASTInterface)parameter).getID() + datatype.getCharacter(i)+"<=>" + datatype.getCharacter(j) + "\t");
					}
					
				}else {
					
					// Asymmetric
					for (int j = 0; j < datatype.getStateCount(); j++) {
						if (i != j) {
							out.append(((BEASTInterface)parameter).getID() + datatype.getCharacter(i)+"=>" + datatype.getCharacter(j) + "\t");
						}
					}
					
				}
				
			}
		
		}else {
			
			for (int i = 0; i < parameterInput.get().getDimension(); i++) {
				out.append(((BEASTInterface)parameter).getID() + "." + (i+1) + "\t");
			}
			
		}
		
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		//parameter.log(sample, out);
		

		
		final RealParameter var = (RealParameter) parameter.getCurrent();
        final int dim = var.getDimension();
        for (int i = 0; i < dim; i++) {
        	double val = var.getValue(i);
        	if (indicatorInput.get() != null && indicatorInput.get().getValue(i) == false) {
        		val = 0.0;
        	}
            out.print(val + "\t");
        }
      
		
	}
	
	@Override
	public void close(PrintStream out) {
		
	}

}





