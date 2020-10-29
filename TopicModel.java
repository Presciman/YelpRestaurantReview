import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Random;

public abstract class TopicModel {

	public String inputFilename;
	public boolean burnedIn = false;

	public void run(int iters, int samples, String filename) throws Exception {
		inputFilename = filename;
		readDocs(filename);
		initialize();
		
		System.out.println("Sampling...");
		
		for (int iter = 1; iter <= iters; iter++) {
			if (iter > (iters-samples)) burnedIn = true;

			System.out.println("Iteration "+iter);
			doSampling(iter);

			// save the output periodically
			/*if (iter % 500 == 0) {// || (iter >= 3000 && iter % 100 == 0)) {
			    System.out.println("Saving output...");
			    writeOutput(filename+iter);
			}*/
		}
		
		// write variable assignments

		writeOutput(filename);

		System.out.println("...done.");
	}
	
	public abstract void initialize();
	
	public abstract void doSampling(int iter);
	
	public abstract void readDocs(String filename) throws Exception;
	
	public abstract void writeOutput(String filename) throws Exception;
}
