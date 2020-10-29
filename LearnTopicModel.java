import java.util.HashMap;

public class LearnTopicModel {

	public static HashMap<String,String> arguments;
	
	public static void main(String[] args) throws Exception {
		arguments = new HashMap<String,String>();
		
		for (int i = 0; i < args.length; i += 2) {
			arguments.put(args[i], args[i+1]);
		}

		String model = arguments.get("-model");
		String filename = arguments.get("-input");
		
		if (model == null) {
			System.out.println("No model specified.");
			return;
		}
		
		if (filename == null) {
			System.out.println("No input file given.");
			return;
		}
		
		TopicModel topicModel = null;

		if (model.equals("flda")) {
			if (!arguments.containsKey("-Z")) {
				System.out.println("Must specify number of topics using -Z");
				return;
			}
			if (!arguments.containsKey("-Y")) {
				System.out.println("Must specify number of aspects using -Y");
				return;
			}
			if (!arguments.containsKey("-K")) {
				System.out.println("Must specify number of dimensions using -K");
				return;
			}
			
			int Z = Integer.parseInt(arguments.get("-Z"));
			int Y = Integer.parseInt(arguments.get("-Y"));
			int K = Integer.parseInt(arguments.get("-K"));
			
			int[] z = new int[K];
			z[0] = Z;
			z[1] = Y;
			for (int i = 2; i < K; i++) z[i] = Y;

			double sigmaA = 1.0;
			if (arguments.containsKey("-sigmaA")) 
				sigmaA = Double.parseDouble(arguments.get("-sigmaA"));
			double sigmaAB = 1.0;
			if (arguments.containsKey("-sigmaAB")) 
				sigmaAB = Double.parseDouble(arguments.get("-sigmaAB"));
			double sigmaW = 0.5;
			if (arguments.containsKey("-sigmaW")) 
				sigmaW = Double.parseDouble(arguments.get("-sigmaW"));
			double sigmaWB = 10.0;
			if (arguments.containsKey("-sigmaWB")) 
				sigmaWB = Double.parseDouble(arguments.get("-sigmaWB"));
			double delta0 = 0.1;
			if (arguments.containsKey("-delta0")) 
				delta0 = Double.parseDouble(arguments.get("-delta0"));
			double delta1 = 0.1;
			if (arguments.containsKey("-delta1")) 
				delta1 = Double.parseDouble(arguments.get("-delta1"));
			double alphaB = -5.0;
			if (arguments.containsKey("-alphaB")) 
				alphaB = Double.parseDouble(arguments.get("-alphaB"));
			double omegaB = -5.0;
			if (arguments.containsKey("-omegaB")) 
				omegaB = Double.parseDouble(arguments.get("-omegaB"));
			double stepSizeADZ = 1e-2;
			if (arguments.containsKey("-stepSizeADZ")) 
				stepSizeADZ = Double.parseDouble(arguments.get("-stepSizeADZ"));
			double stepSizeAZ = stepSizeADZ / 100.0;
			if (arguments.containsKey("-stepSizeAZ")) 
				stepSizeAZ = Double.parseDouble(arguments.get("-stepSizeAZ"));
			double stepSizeAB = stepSizeADZ / 100.0;
			if (arguments.containsKey("-stepSizeAB")) 
				stepSizeAB = Double.parseDouble(arguments.get("-stepSizeAB"));
			double stepSizeW = 1e-3;
			if (arguments.containsKey("-stepSizeW")) 
				stepSizeW = Double.parseDouble(arguments.get("-stepSizeW"));
			double stepSizeWB = stepSizeW / 100.0;
			if (arguments.containsKey("-stepSizeWB")) 
				stepSizeWB = Double.parseDouble(arguments.get("-stepSizeWB"));
			double stepSizeB = 1e-3;
			if (arguments.containsKey("-stepSizeB")) 
				stepSizeB = Double.parseDouble(arguments.get("-stepSizeB"));

			int likelihoodFreq = 100;
			if (arguments.containsKey("-likelihoodFreq")) 
				likelihoodFreq = Integer.parseInt(arguments.get("-likelihoodFreq"));
			if (likelihoodFreq == 0) likelihoodFreq = 1;
			else if (likelihoodFreq == -1) likelihoodFreq = Integer.MAX_VALUE; 
			else if (likelihoodFreq < -1) {
				System.out.println("Invalid value for likelihoodFreq; must be positive");
				return;	
			}
			int blockFreq = 0;
			if (arguments.containsKey("-blockFreq")) 
				blockFreq = Integer.parseInt(arguments.get("-blockFreq"));
			if (blockFreq == 0) blockFreq = 1;
			else if (blockFreq == -1) blockFreq = Integer.MAX_VALUE; 
			else if (blockFreq < -1) {
				System.out.println("Invalid value for blockFreq; must be positive");
				return;	
			}

			String priorPrefix = "";
			if (arguments.containsKey("-priorPrefix")) 
				priorPrefix = arguments.get("-priorPrefix");

			topicModel = new FLDA(z, sigmaA, sigmaAB, sigmaW, sigmaWB, stepSizeADZ, stepSizeAZ, stepSizeAB, stepSizeW, stepSizeWB, stepSizeB, delta0, delta1, alphaB, omegaB, likelihoodFreq, blockFreq, priorPrefix); 
		}
		else {
			System.out.println("Invalid model specification. Options: flda");
			return;
		}
		
		int iters = 2000;
		if (arguments.containsKey("-iters")) 
			iters = Integer.parseInt(arguments.get("-iters"));
		int samples = 100;
		if (arguments.containsKey("-samples")) 
			samples = Integer.parseInt(arguments.get("-samples"));
		
		topicModel.run(iters, samples, filename);
	}

}
