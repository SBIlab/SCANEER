package scripts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import covariance.datacontainers.Alignment;
import covariance.parsers.PfamParser;
import covariance.utils.ConfigReader;

public class DumpPfamNames
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				ConfigReader.getOutDataDir() + File.separator + "pfamNames.txt"
				)));
		
		writer.write("pfamId\tnumSequences\tnumColumns\n");
		
		PfamParser parser = new PfamParser();
		
		while ( ! parser.isFinished())
		{
			Alignment a = parser.getNextAlignment();
			
			if ( a != null)
			{
				writer.write(a.getAligmentID() + "\t" + a.getNumSequencesInAlignment() + "\t" + 
						a.getNumColumnsInAlignment() + "\n");
				System.out.println(a.getAligmentID());
				writer.flush();
			}
		}
		
		writer.flush();  writer.close();
	}
}
