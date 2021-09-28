package covariance.utils;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

public class ProcessWrapper
{
	Process process;
	InputStream in;
	InputStream inError;
	
	public ProcessWrapper( String[] cmdArgs ) throws Exception
	{
		Runtime r = Runtime.getRuntime();
		Process p = r.exec(cmdArgs);
		
		BufferedReader br = new BufferedReader (new InputStreamReader(p.getInputStream ()));
		
		String s;
		
		while ((s = br.readLine ())!= null)
		{
    		//System.out.println (s);
		}
				
		p.waitFor();
		p.destroy();
	}
}
