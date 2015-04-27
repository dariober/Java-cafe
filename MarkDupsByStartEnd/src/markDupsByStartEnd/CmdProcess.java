package markDupsByStartEnd;

import java.io.InputStream;
import java.io.OutputStream;

// TODO: Allow to store the string with the executed command.
public class CmdProcess extends Process{

	private static String cmd= "";

	public static String getCmd() {
		return cmd;
	}

	public static void setCmd(String cmd) {
		CmdProcess.cmd = cmd;
	}

	@Override
	public void destroy() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int exitValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public InputStream getErrorStream() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public InputStream getInputStream() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public OutputStream getOutputStream() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int waitFor() throws InterruptedException {
		// TODO Auto-generated method stub
		return 0;
	}
	
}
